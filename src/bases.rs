use std::string::String;

#[derive(Debug, PartialEq)]
pub struct Bases(String);

impl Bases {
    #[allow(dead_code)]
    pub fn new(bases: String) -> Bases {
        Bases(bases)
    }

    pub fn from_str(bases: &str) -> Bases {
        Bases(bases.to_string())
    }
}
