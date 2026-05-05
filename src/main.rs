use std::fs::File;
use std::io::Read;
use swl::parse;
use std::env;

fn read_file(path: &str) -> Result<String, std::io::Error> {
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

fn main() {
    let argv: Vec<String> = env::args().collect();
    let code = read_file(&argv[1]).ok();
    if let Some(code) = code {
        let result = parse(code.as_bytes());
        println!("{}", result);
    } else {
        println!("Could not code from file");
    }
}
