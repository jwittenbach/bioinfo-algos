use std::path::PathBuf;
use std::fs::read_to_string;
use std::convert::TryFrom;
use phf::phf_map;

const BASE_MAP: phf::Map<char, u8> = phf_map! {
    'a' => 0,
    'A' => 0,
    'c' => 1,
    'C' => 1,
    'g' => 2,
    'G' => 2,
    't' => 3,
    'T' => 3
};

const INDEX_MAP: [char; 4] = ['A', 'C', 'G', 'T'];

pub fn from_str(s: &str) -> Vec<u8> {
    s.chars().map(|c| BASE_MAP.get(&c).unwrap().clone()).collect()
}

pub fn from_file(path: &PathBuf) -> Vec<u8> {
    from_str(&read_to_string(&path).unwrap().trim())
}

pub fn to_str(seq: &[u8]) -> String {
    seq.iter().map(|&i| INDEX_MAP[i as usize]).collect()
}

pub fn revc(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|x| 3 - x).collect()
}

pub fn seq_to_num(pattern: &[u8]) -> usize {
    pattern.iter().fold(0, |acc, &i| (4 * acc) + usize::try_from(i).unwrap())
}

pub fn num_to_seq(num: usize, k: u32) -> Vec<u8> {
    let mut n = num;
    let mut v = Vec::new();
    for _ in 0..k {
        v.push((n % 4) as u8);
        n = n / 4;
    }
    v.into_iter().rev().collect()
}