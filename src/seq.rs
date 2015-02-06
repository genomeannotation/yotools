use std::string::String;

#[derive(Debug, PartialEq)]
pub struct Bases(String);

impl Bases {
    pub fn new(bases: String) -> Bases {
        Bases(bases)
    }
}
