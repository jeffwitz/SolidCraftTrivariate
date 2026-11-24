//! JSON interchange for the offline demonstration CLI. Number tokens retain
//! their spelling so exact integer coefficients never pass through binary64.

use std::collections::BTreeMap;
use std::fmt::{self, Write};

#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Number(String),
    String(String),
    Array(Vec<Json>),
    Object(BTreeMap<String, Json>),
}

pub fn object<const N: usize>(pairs: [(&str, Json); N]) -> Json {
    Json::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

impl Json {
    pub fn text(&self) -> Result<&str, &'static str> {
        match self {
            Self::String(s) => Ok(s),
            _ => Err("expected a JSON string"),
        }
    }
    pub fn items(&self) -> Result<&[Self], &'static str> {
        match self {
            Self::Array(a) => Ok(a),
            _ => Err("expected a JSON array"),
        }
    }
    pub fn field(&self, key: &str) -> Result<&Self, &'static str> {
        match self {
            Self::Object(o) => o.get(key).ok_or("missing object field"),
            _ => Err("expected a JSON object"),
        }
    }
    pub fn keys(&self, keys: &[&str]) -> Result<(), &'static str> {
        match self {
            Self::Object(o) if o.len() == keys.len() && keys.iter().all(|k| o.contains_key(*k)) => {
                Ok(())
            }
            _ => Err("unexpected object fields"),
        }
    }
    pub fn number<T: std::str::FromStr>(&self) -> Result<T, &'static str> {
        match self {
            Self::Number(s) => s.parse().map_err(|_| "number outside expected type"),
            _ => Err("expected a JSON number"),
        }
    }
}

fn quoted(s: &str, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    f.write_char('"')?;
    for c in s.chars() {
        match c {
            '"' => f.write_str("\\\"")?,
            '\\' => f.write_str("\\\\")?,
            '\n' => f.write_str("\\n")?,
            '\r' => f.write_str("\\r")?,
            '\t' => f.write_str("\\t")?,
            c if c < ' ' => write!(f, "\\u{:04x}", c as u32)?,
            c => f.write_char(c)?,
        }
    }
    f.write_char('"')
}

impl fmt::Display for Json {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Null => f.write_str("null"),
            Self::Bool(b) => write!(f, "{b}"),
            Self::Number(n) => f.write_str(n),
            Self::String(s) => quoted(s, f),
            Self::Array(a) => {
                f.write_char('[')?;
                for (i, v) in a.iter().enumerate() {
                    if i > 0 {
                        f.write_char(',')?;
                    }
                    write!(f, "{v}")?;
                }
                f.write_char(']')
            }
            Self::Object(o) => {
                f.write_char('{')?;
                for (i, (k, v)) in o.iter().enumerate() {
                    if i > 0 {
                        f.write_char(',')?;
                    }
                    quoted(k, f)?;
                    write!(f, ":{v}")?;
                }
                f.write_char('}')
            }
        }
    }
}

struct Parser<'a>(&'a str);

impl Parser<'_> {
    fn ws(&mut self) {
        self.0 = self.0.trim_start_matches([' ', '\t', '\r', '\n']);
    }
    fn eat(&mut self, s: &str) -> bool {
        if let Some(rest) = self.0.strip_prefix(s) {
            self.0 = rest;
            true
        } else {
            false
        }
    }
    fn next(&mut self) -> Result<char, &'static str> {
        let c = self.0.chars().next().ok_or("unexpected end of JSON")?;
        self.0 = &self.0[c.len_utf8()..];
        Ok(c)
    }
    fn hex4(&mut self) -> Result<u32, &'static str> {
        let mut value = 0;
        for _ in 0..4 {
            value = 16 * value + self.next()?.to_digit(16).ok_or("invalid Unicode escape")?;
        }
        Ok(value)
    }
    fn string(&mut self) -> Result<String, &'static str> {
        if !self.eat("\"") {
            return Err("expected string");
        }
        let mut result = String::new();
        loop {
            let c = self.next()?;
            match c {
                '"' => return Ok(result),
                c if c < ' ' => return Err("unescaped control character"),
                '\\' => result.push(match self.next()? {
                    c @ ('"' | '\\' | '/') => c,
                    'b' => '\u{8}',
                    'f' => '\u{c}',
                    'n' => '\n',
                    'r' => '\r',
                    't' => '\t',
                    'u' => {
                        let mut code = self.hex4()?;
                        if (0xd800..=0xdbff).contains(&code) {
                            if !self.eat("\\u") {
                                return Err("missing low surrogate");
                            }
                            let low = self.hex4()?;
                            if !(0xdc00..=0xdfff).contains(&low) {
                                return Err("invalid surrogate pair");
                            }
                            code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                        }
                        char::from_u32(code).ok_or("invalid Unicode scalar")?
                    }
                    _ => return Err("invalid string escape"),
                }),
                c => result.push(c),
            }
        }
    }
    fn digits(&mut self) -> Result<(), &'static str> {
        let count = self.0.bytes().take_while(u8::is_ascii_digit).count();
        if count == 0 {
            return Err("expected digits");
        }
        self.0 = &self.0[count..];
        Ok(())
    }
    fn numeric(&mut self) -> Result<Json, &'static str> {
        let start = self.0;
        self.eat("-");
        if !self.eat("0") {
            self.digits()?;
        }
        if self.eat(".") {
            self.digits()?;
        }
        if self.eat("e") || self.eat("E") {
            if !self.eat("-") {
                self.eat("+");
            }
            self.digits()?;
        }
        Ok(Json::Number(start[..start.len() - self.0.len()].into()))
    }
    fn value(&mut self, depth: usize) -> Result<Json, &'static str> {
        if depth > 128 {
            return Err("JSON nesting exceeds limit");
        }
        self.ws();
        if self.eat("{") {
            let mut map = BTreeMap::new();
            self.ws();
            if self.eat("}") {
                return Ok(Json::Object(map));
            }
            loop {
                self.ws();
                let key = self.string()?;
                self.ws();
                if !self.eat(":") {
                    return Err("expected colon");
                }
                let val = self.value(depth + 1)?;
                if map.insert(key, val).is_some() {
                    return Err("duplicate JSON object key");
                }
                self.ws();
                if self.eat("}") {
                    break;
                }
                if !self.eat(",") {
                    return Err("expected comma");
                }
            }
            Ok(Json::Object(map))
        } else if self.eat("[") {
            let mut arr = vec![];
            self.ws();
            if self.eat("]") {
                return Ok(Json::Array(arr));
            }
            loop {
                arr.push(self.value(depth + 1)?);
                self.ws();
                if self.eat("]") {
                    break;
                }
                if !self.eat(",") {
                    return Err("expected comma");
                }
            }
            Ok(Json::Array(arr))
        } else if self.0.starts_with('"') {
            Ok(Json::String(self.string()?))
        } else if self.eat("true") {
            Ok(Json::Bool(true))
        } else if self.eat("false") {
            Ok(Json::Bool(false))
        } else if self.eat("null") {
            Ok(Json::Null)
        } else {
            self.numeric()
        }
    }
}

pub fn parse(source: &str) -> Result<Json, &'static str> {
    let mut p = Parser(source);
    let value = p.value(0)?;
    p.ws();
    if !p.0.is_empty() {
        return Err("trailing JSON content");
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_numbers_unicode_and_invalid_tokens() {
        let value = parse(r#"{"a":[1234567891234567891,-1.25e+2,"\uD834\uDD1E\n"]}"#).unwrap();
        assert_eq!(
            value.field("a").unwrap().items().unwrap()[0]
                .number::<i64>()
                .unwrap(),
            1234567891234567891
        );
        assert_eq!(parse(&value.to_string()), Ok(value));
        for invalid in [
            "1.",
            "01",
            "[2,]",
            "{",
            "2e-",
            "{\"a\":1,\"a\":2}",
            "\"\\uD800\"",
            "NaN",
        ] {
            assert!(parse(invalid).is_err(), "{invalid}");
        }
    }
}
