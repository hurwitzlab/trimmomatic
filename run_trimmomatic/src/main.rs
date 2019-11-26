extern crate run_trimmomatic;
use std::process;

fn main() {
    let config = match run_trimmomatic::get_args() {
        Ok(c) => c,
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };

    if let Err(e) = run_trimmomatic::run(config) {
        println!("Error: {}", e);
        process::exit(1);
    }
}
