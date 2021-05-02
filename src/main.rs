mod seqs {
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

    pub fn to_str(seq: &[u8]) -> String {
        seq.iter().map(|&i| INDEX_MAP[i as usize]).collect()
    }
}

//------------------------------
fn pattern_count(text: &[u8], pattern: &[u8]) -> usize {
    text.windows(pattern.len()).fold(0, |acc, kmer| acc + (if kmer == pattern {1} else {0}))
}

use std::cmp::max;
use std::collections::HashSet;

fn frequent_words(text: &[u8], k: usize) -> HashSet<&[u8]> {
    let count_kmer_pairs: Vec<_> = text.windows(k).map(|kmer| (kmer, pattern_count(text, kmer))).collect();
    let max_count = count_kmer_pairs.iter().fold(0, |acc, &(_, count)| max(acc, count));
    count_kmer_pairs.iter().filter(|(_, count)| max_count == *count).map(|&(kmer, _)| kmer).collect()
}

//------------------------------

mod cli {
    use structopt::StructOpt;

    #[derive(StructOpt, Debug)]
    pub enum Command {

        CountPattern {
            #[structopt(short)]
            text: String,
            #[structopt(short)]
            pattern: String
        },

        FrequentWords {
            #[structopt(short)]
            text: String,
            #[structopt(short)]
            k: usize 
        }
    }

}

use structopt::StructOpt;

fn main() {
    match cli::Command::from_args() {

        cli::Command::CountPattern{text, pattern} => println!("{:?}",
            pattern_count(&seqs::from_str(&text[..]), &seqs::from_str(&pattern[..]))
        ),

        cli::Command::FrequentWords{text, k} => println!("{:?}",
            frequent_words(&seqs::from_str(&text[..]), k).iter().map(|seq| seqs::to_str(&seq)).collect::<Vec<_>>()
        )
    }
}
