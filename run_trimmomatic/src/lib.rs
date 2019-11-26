extern crate clap;

use clap::{App, Arg};
use regex::Regex;
use std::collections::HashMap;
use std::error::Error;
use std::process::{Command, Stdio};
use std::{
    env,
    fs::{self, DirBuilder},
    io::Write,
    path::{Path, PathBuf},
};

#[derive(Debug)]
pub struct Config {
    query: Vec<String>,
    trimmomatic_jar: String,
    out_dir: PathBuf,
    phred_base: Option<PhredBase>,
    num_concurrent_jobs: Option<u32>,
    num_halt: Option<u32>,
    num_threads: Option<u32>,
}

#[derive(Debug)]
struct SplitPath {
    stem: String,
    ext: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Hash)]
enum ReadDirection {
    Forward,
    Reverse,
}

#[derive(Debug)]
enum PhredBase {
    Phred33,
    Phred64,
}

type MyResult<T> = Result<T, Box<dyn Error>>;
type ReadPair = HashMap<ReadDirection, String>;
type ReadPairLookup = HashMap<String, ReadPair>;
type SingleReads = Vec<String>;

// --------------------------------------------------
pub fn get_args() -> MyResult<Config> {
    let matches = App::new("run_trimmomatic")
        .version("0.1.0")
        .author("Ken Youens-Clark <kyclark@email.arizona.edu>")
        .about("Runs Pear")
        .arg(
            Arg::with_name("query")
                .short("q")
                .long("query")
                .value_name("FILE_OR_DIR")
                .help("File input or directory")
                .required(true)
                .min_values(1),
        )
        .arg(
            Arg::with_name("trimmomatic_jar")
                .short("t")
                .long("trimmomatic")
                .value_name("STR")
                .required(true)
                .help("Location of Trimmomatic jar file"),
        )
        .arg(
            Arg::with_name("out_dir")
                .short("o")
                .long("out_dir")
                .value_name("DIR")
                .help("Output directory"),
        )
        .arg(
            Arg::with_name("phred33")
                .long("phred33")
                .help("Phread33 base")
                .conflicts_with("phred64"),
        )
        .arg(
            Arg::with_name("phred64")
                .long("phred64")
                .help("Phread64 base")
                .conflicts_with("phred33"),
        )
        .arg(
            Arg::with_name("threads")
                .short("j")
                .long("threads")
                .value_name("INT")
                .help("Number of threads to use"),
        )
        .arg(
            Arg::with_name("num_concurrent_jobs")
                .short("J")
                .long("num_concurrent_jobs")
                .value_name("INT")
                .default_value("8")
                .help("Number of concurrent jobs for parallel"),
        )
        .arg(
            Arg::with_name("num_halt")
                .short("H")
                .long("num_halt")
                .value_name("INT")
                .default_value("1")
                .help("Halt after this many failing jobs"),
        )
        .get_matches();

    let out_dir = match matches.value_of("out_dir") {
        Some(x) => PathBuf::from(x),
        _ => {
            let cwd = env::current_dir()?;
            cwd.join(PathBuf::from("trimmomatic-out"))
        }
    };

    let num_concurrent_jobs = matches
        .value_of("num_concurrent_jobs")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let num_halt = matches
        .value_of("num_halt")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let num_threads = matches
        .value_of("threads")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let phred_base = if matches.is_present("phred33") {
        Some(PhredBase::Phred33)
    } else if matches.is_present("phred64") {
        Some(PhredBase::Phred64)
    } else {
        None
    };

    Ok(Config {
        query: matches.values_of_lossy("query").unwrap(),
        trimmomatic_jar: matches
            .value_of_lossy("trimmomatic_jar")
            .unwrap()
            .to_string(),
        out_dir,
        phred_base,
        num_concurrent_jobs,
        num_halt,
        num_threads,
    })
}

// --------------------------------------------------
pub fn run(config: Config) -> MyResult<()> {
    println!("config {:?}", config);
    let files = find_files(&config.query)?;

    if files.is_empty() {
        let msg = format!("No input files from query \"{:?}\"", &config.query);
        return Err(From::from(msg));
    }

    let (pairs, singles) = classify(&files)?;

    println!(
        "Processing {} pair, {} single.",
        pairs.keys().len(),
        singles.len()
    );

    let jobs = make_jobs(&config, pairs, singles)?;

    println!("{:?}", jobs);
    run_jobs(
        &jobs,
        "Running pear",
        config.num_concurrent_jobs.unwrap_or(8),
        config.num_halt.unwrap_or(1),
    )?;

    println!("Done, see output in \"{}\"", &config.out_dir.display());

    Ok(())
}

// --------------------------------------------------
fn find_files(paths: &[String]) -> Result<Vec<String>, Box<dyn Error>> {
    let mut files = vec![];
    for path in paths {
        let meta = fs::metadata(path)?;
        if meta.is_file() {
            files.push(path.to_owned());
        } else {
            for entry in fs::read_dir(path)? {
                let entry = entry?;
                let meta = entry.metadata()?;
                if meta.is_file() {
                    files.push(entry.path().display().to_string());
                }
            }
        };
    }

    if files.is_empty() {
        return Err(From::from("No input files"));
    }

    Ok(files)
}

// --------------------------------------------------
fn classify(
    paths: &[String],
) -> Result<(ReadPairLookup, SingleReads), Box<dyn Error>> {
    let paths = paths.iter().map(Path::new);
    let mut exts: Vec<String> =
        paths.clone().map(get_extension).filter_map(|x| x).collect();
    exts.dedup();

    let dots = Regex::new(r"\.").unwrap();
    let exts: Vec<String> = exts
        .into_iter()
        .map(|x| dots.replace(&x, r"\.").to_string())
        .collect();

    let pattern = format!(r"(.+)[_-][Rr]?([12])?\.(?:{})$", exts.join("|"));
    let re = Regex::new(&pattern).unwrap();
    let mut pairs: ReadPairLookup = HashMap::new();
    let mut singles: Vec<String> = vec![];

    for path in paths.map(Path::new) {
        let path_str = path.to_str().expect("Convert path");

        if let Some(file_name) = path.file_name() {
            let basename = file_name.to_string_lossy();
            if let Some(cap) = re.captures(&basename) {
                let sample_name = &cap[1];
                let direction = if &cap[2] == "1" {
                    ReadDirection::Forward
                } else {
                    ReadDirection::Reverse
                };

                if !pairs.contains_key(sample_name) {
                    let mut pair: ReadPair = HashMap::new();
                    pair.insert(direction, path_str.to_string());
                    pairs.insert(sample_name.to_string(), pair);
                } else if let Some(pair) = pairs.get_mut(sample_name) {
                    pair.insert(direction, path_str.to_string());
                }
            } else {
                singles.push(path_str.to_string());
            }
        }
    }

    let bad: Vec<String> = pairs
        .iter()
        .filter_map(|(k, v)| {
            if !v.contains_key(&ReadDirection::Forward)
                || !v.contains_key(&ReadDirection::Reverse)
            {
                Some(k.to_string())
            } else {
                None
            }
        })
        .collect();

    // Push unpaired samples to the singles
    for key in bad {
        pairs.remove(&key);
        singles.push(key);
    }

    Ok((pairs, singles))
}

// --------------------------------------------------
/// Returns the extension plus optional ".gz"
fn get_extension(path: &Path) -> Option<String> {
    if let Some(basename) = path.file_name() {
        let res = split_filename(basename.to_string_lossy().to_string());
        return res.ext;
    }
    None
}

// --------------------------------------------------
fn split_filename(filename: String) -> SplitPath {
    let re = Regex::new(r"([^/]+?)\.([^.]+(?:\.gz)?)$").unwrap();
    if let Some(cap) = re.captures(&filename) {
        SplitPath {
            stem: cap[1].to_string(),
            ext: Some(cap[2].to_string()),
        }
    } else {
        SplitPath {
            stem: filename,
            ext: None,
        }
    }
}

// --------------------------------------------------
fn make_jobs(
    config: &Config,
    pairs: ReadPairLookup,
    singles: SingleReads,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut args: Vec<String> = vec![];

    if let Some(phred_base) = &config.phred_base {
        args.push(format!(
            "{}",
            match phred_base {
                PhredBase::Phred33 => "-phred33".to_string(),
                _ => "-phred64".to_string(),
            }
        ));
    }

    println!("pairs {:?}", pairs);
    println!("singles {:?}", singles);

    let mut jobs: Vec<String> = vec![];
    for (i, (sample, val)) in pairs.iter().enumerate() {
        println!("{:3}: Pair {}", i + 1, sample);

        if let (Some(fwd), Some(rev)) = (
            val.get(&ReadDirection::Forward),
            val.get(&ReadDirection::Reverse),
        ) {
            let out_dir = &config.out_dir.join(sample);
            if !out_dir.is_dir() {
                DirBuilder::new().recursive(true).create(&out_dir)?;
            }

            let SplitPath {
                stem: fwd_stem,
                ext: fwd_ext,
            } = split_filename(fwd.to_string());

            let fwd_ext = &fwd_ext.unwrap_or("".to_string());

            let SplitPath {
                stem: rev_stem,
                ext: rev_ext,
            } = split_filename(rev.to_string());

            let rev_ext = &rev_ext.unwrap_or("".to_string());

            let fwd_out = &out_dir.join(format!("{}.{}", fwd_stem, fwd_ext));
            let fwd_out_u =
                &out_dir.join(format!("{}.unpaired.{}", fwd_stem, &fwd_ext));
            let rev_out = &out_dir.join(format!("{}.{}", rev_stem, &rev_ext));
            let rev_out_u =
                &out_dir.join(format!("{}.unpaired.{}", rev_stem, &rev_ext));
            let log_file = out_dir.join(format!("{}.log", sample));

            jobs.push(format!(
                "java -jar {} PE -trimlog {} {} {} \
                 {} {} {} {} {} SLIDINGWINDOW:40:15",
                &config.trimmomatic_jar,
                log_file.display(),
                fwd,
                rev,
                fwd_out.display(),
                fwd_out_u.display(),
                rev_out.display(),
                rev_out_u.display(),
                args.join(" "),
            ));
        }
    }

    for (i, file) in singles.iter().enumerate() {
        println!("{:3}: Single {}", i + 1, file);

        let path = Path::new(file);
        let basename = path.file_name().expect("basename");
        let SplitPath { stem, ext: _ } =
            split_filename(basename.to_string_lossy().to_string());
        let out_dir = &config.out_dir.join(&stem);
        if !out_dir.is_dir() {
            DirBuilder::new().recursive(true).create(&out_dir)?;
        }

        let out_file = out_dir.join(&basename);
        if let Some(log_file) =
            out_dir.join(format!("{}.log", stem)).as_path().to_str()
        {
            args.push(format!("-trimlog {}", log_file));
        }

        jobs.push(format!(
            "java -jar {} SE {} {} {} SLIDINGWINDOW:40:15",
            &config.trimmomatic_jar,
            args.join(" "),
            file,
            out_file.display(),
        ));
    }

    Ok(jobs)
}

// --------------------------------------------------
fn run_jobs(
    jobs: &[String],
    msg: &str,
    num_concurrent_jobs: u32,
    num_halt: u32,
) -> MyResult<()> {
    let num_jobs = jobs.len();

    if num_jobs > 0 {
        println!(
            "{} (# {} job{} @ {})",
            msg,
            num_jobs,
            if num_jobs == 1 { "" } else { "s" },
            num_concurrent_jobs,
        );

        let mut args: Vec<String> =
            vec!["-j".to_string(), num_concurrent_jobs.to_string()];

        if num_halt > 0 {
            args.push("--halt".to_string());
            args.push(format!("soon,fail={}", num_halt.to_string()));
        }

        let mut process = Command::new("parallel")
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .spawn()?;

        {
            let stdin = process.stdin.as_mut().expect("Failed to open stdin");
            stdin
                .write_all(jobs.join("\n").as_bytes())
                .expect("Failed to write to stdin");
        }

        let result = process.wait()?;
        if !result.success() {
            return Err(From::from("Failed to run jobs in parallel"));
        }
    }

    Ok(())
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_extension() {
        assert_eq!(
            get_extension(Path::new("foo.fna")),
            Some("fna".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fasta.gz")),
            Some("fasta.gz".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fa.gz")),
            Some("fa.gz".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fasta")),
            Some("fasta".to_string())
        );

        assert_eq!(get_extension(Path::new("foo.fq")), Some("fq".to_string()));

        assert_eq!(get_extension(Path::new("foo")), None);
    }

    #[test]
    fn test_classify() {
        let res = classify(&["ERR1711926.fastq.gz".to_string()]);
        assert!(res.is_ok());

        if let Ok((pairs, singles)) = res {
            assert_eq!(pairs.len(), 0);
            assert_eq!(singles.len(), 1);
        }

        let res = classify(&[
            "/foo/bar/ERR1711926_1.fastq.gz".to_string(),
            "/foo/bar/ERR1711926_2.fastq.gz".to_string(),
            "/foo/bar/ERR1711927-R1.fastq.gz".to_string(),
            "/foo/bar/ERR1711927_R2.fastq.gz".to_string(),
            "/foo/bar/ERR1711928.fastq.gz".to_string(),
            "/foo/bar/ERR1711929_1.fastq.gz".to_string(),
        ]);
        assert!(res.is_ok());

        if let Ok((pairs, singles)) = res {
            assert_eq!(pairs.len(), 2);
            assert_eq!(singles.len(), 2);

            assert!(pairs.contains_key("ERR1711926"));
            assert!(pairs.contains_key("ERR1711927"));

            //assert!(!singles.contains_key("ERR1711928"));
            //assert!(!singles.contains_key("ERR1711929"));

            if let Some(val) = pairs.get("ERR1711926") {
                assert!(val.contains_key(&ReadDirection::Forward));
                assert!(val.contains_key(&ReadDirection::Reverse));

                if let Some(fwd) = val.get(&ReadDirection::Forward) {
                    assert_eq!(fwd, &"/foo/bar/ERR1711926_1.fastq.gz");
                }
                if let Some(rev) = val.get(&ReadDirection::Reverse) {
                    assert_eq!(rev, &"/foo/bar/ERR1711926_2.fastq.gz");
                }
            }

            if let Some(val) = pairs.get("ERR1711927") {
                assert!(val.contains_key(&ReadDirection::Forward));
                assert!(val.contains_key(&ReadDirection::Reverse));

                if let Some(fwd) = val.get(&ReadDirection::Forward) {
                    assert_eq!(fwd, &"/foo/bar/ERR1711927-R1.fastq.gz");
                }
                if let Some(rev) = val.get(&ReadDirection::Reverse) {
                    assert_eq!(rev, &"/foo/bar/ERR1711927_R2.fastq.gz");
                }
            }
        }
    }

    #[test]
    fn test_split_filename() {
        let res = split_filename("foo.fa".to_string());
        assert_eq!(res.stem, "foo".to_string());
        assert_eq!(res.ext, Some("fa".to_string()));

        let res = split_filename("foo.fa.gz".to_string());
        assert_eq!(res.stem, "foo".to_string());
        assert_eq!(res.ext, Some("fa.gz".to_string()));

        let res = split_filename("/foo/bar/baz.fa.gz".to_string());
        assert_eq!(res.stem, "baz".to_string());
        assert_eq!(res.ext, Some("fa.gz".to_string()));
    }
}
