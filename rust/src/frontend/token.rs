
#[derive(partialEq, Debug, Clone)]
pub enum Token {
    EOF,
    Identifier(String),
    StringLiteral(String),
    Assign,
    Pipe,
    Colon,
    Comma,
    Lambda,
    Arrow,
    LParen,
    RParen,
    LBrace,
    RBrace,
    Import,
}

#[derive(Clone, Copy, PartialEq, Debug)]
#[repr(C)]
pub struct Tokens<'a> {
    pub rep: &'a [Token],
    pub start: usize,
    pub end: usize
}

