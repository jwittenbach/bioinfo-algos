mod seqs {
    use std::path::PathBuf;
    use std::fs::read_to_string;
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

}

//------------------------------
fn pattern_count(text: &[u8], pattern: &[u8]) -> usize {
    text.windows(pattern.len()).fold(0, |acc, kmer| acc + (if kmer == pattern {1} else {0}))
}

use std::cmp::max;
use std::convert::TryFrom;
use std::collections::{HashMap, HashSet};

fn seq_to_num(pattern: &[u8]) -> usize {
    pattern.iter().fold(0, |acc, &i| (4 * acc) + usize::try_from(i).unwrap())
}

fn num_to_seq(num: usize, k: u32) -> Vec<u8> {
    let mut n = num;
    let mut v = Vec::new();
    for _ in 0..k {
        v.push((n % 4) as u8);
        n = n / 4;
    }
    v.into_iter().rev().collect()
}

fn count_array(text: &[u8], k: u32) -> Vec<u32> {
    let mut counts = vec![0; usize::checked_pow(4, k).unwrap()];
    text.windows(usize::try_from(k).unwrap()).for_each(|kmer| {
        let el = counts.get_mut(seq_to_num(kmer)).unwrap();
        *el += 1;
    });
    counts 
}

fn count_map(text: &[u8], k: u32) -> HashMap<usize, u32> {
    let mut counts = HashMap::new();
    text.windows(usize::try_from(k).unwrap()).for_each(|kmer| {
        let count = counts.entry(seq_to_num(kmer)).or_insert(0);
        *count += 1;
    });
    counts
}

fn frequent_words(text: &[u8], k: u32) -> HashSet<Vec<u8>> {
    let count_kmer_pairs: Vec<_> = text
        .windows(usize::try_from(k).unwrap())
        .map(|kmer| (kmer, pattern_count(text, kmer)))
        .collect();
    let max_count = count_kmer_pairs
        .iter()
        .fold(0, |acc, &(_, count)| max(acc, count));
    count_kmer_pairs
        .iter()
        .filter(|(_, count)| max_count == *count).map(|&(kmer, _)| kmer.iter().cloned().collect())
        .collect()
}

fn frequent_fast(text: &[u8], k: u32) -> HashSet<Vec<u8>> {
    let counts = count_array(&text, k);
    let max_count = counts.iter().reduce(|acc, count| max(acc, count)).unwrap();
    counts
        .iter()
        .enumerate()
        .filter_map(|(idx, count)|
            if count == max_count {Some(num_to_seq(idx, k))} else {None})
        .collect() 
}

fn frequent_faster(text: &[u8], k: u32) -> HashSet<Vec<u8>> {
    let counts = count_map(&text, k);
    let max_count = counts.iter().fold(0, |acc, (_, &count)| max(acc, count));
    counts.iter()
        .filter_map(|(&idx, &count)|
            if count == max_count {Some(num_to_seq(idx, k))} else {None})
        .collect()
}

fn frequent_sort(text: &[u8], k: u32) -> HashSet<Vec<u8>> {
    let mut kmer_ids = text
        .windows(usize::try_from(k).unwrap())
        .map(|kmer| seq_to_num(kmer))
        .collect::<Vec<_>>();
    kmer_ids.sort();
    let mut counts = vec![1; text.len()];
    kmer_ids
        .windows(2)
        .enumerate()
        .for_each(|(idx, ids)|
            if &ids[0] == &ids[1] {
                counts[idx + 1] = counts[idx] + 1;
            });
    let max_count = counts.iter().reduce(|acc, count| max(acc, count)).unwrap();
    kmer_ids
        .iter()
        .zip(counts.iter())
        .filter_map(|(kmer_id, count)| if count == max_count {Some(num_to_seq(*kmer_id, k))} else {None})
        .collect()
} 

fn pattern_match(text: &[u8], pattern: &[u8]) -> Vec<u32> {
    text
        .windows(pattern.len())
        .enumerate()
        .filter_map(|(idx, kmer)| if kmer == pattern {Some(u32::try_from(idx).unwrap())} else {None})
        .collect()
}

fn clumps(text: &[u8], k: u32, t: u32, l: u32) -> HashSet<Vec<u8>> {
    let k_ = usize::try_from(k).unwrap();
    let l_ = usize::try_from(l).unwrap();

    // initialize count array + clumps with counts from first window
    let mut clump = vec![false; usize::checked_pow(4, k).unwrap()];
    let mut counts = count_array(&text[..l_], k);
    counts
        .iter()
        .enumerate()
        .for_each(|(idx, &count)| {
            if count >= t {
                clump[idx] = true;
            }
        });

    let mut update_clumps = |init_kmer, final_kmer| {
        counts[seq_to_num(init_kmer)] -= 1;
        let final_ind = seq_to_num(final_kmer);
        counts[final_ind] += 1;
        if counts[final_ind] >= t {
            clump[final_ind] = true;
        }
    };

    text
        .windows(l_)
        .for_each(|window| {
            update_clumps(&window[..k_], &window[l_ - k_..])
        });

    clump
        .iter()
        .enumerate()
        .filter_map(|(idx, &c)| if c {Some(num_to_seq(idx, k))} else {None})
        .collect()

}

fn clumps_hash(text: &[u8], k: u32, t: u32, l: u32) -> HashSet<Vec<u8>> {
    let k_ = usize::try_from(k).unwrap();
    let l_ = usize::try_from(l).unwrap();

    let mut clump: HashSet<usize> = HashSet::new();
    let mut counts = count_map(&text[..l_], k);

    let mut update_clumps = |init_kmer, final_kmer| {
        *counts.get_mut(&seq_to_num(init_kmer)).unwrap() -= 1;
        let final_ind = seq_to_num(final_kmer);
        let final_count = counts.entry(final_ind).or_insert(0);
        *final_count += 1;
        if *final_count >= t {
            clump.insert(final_ind);
        }
    };

    text
        .windows(l_)
        .for_each(|window| {
            update_clumps(&window[..k_], &window[l_ - k_..])
        });

    clump
        .iter()
        .map(|&idx| num_to_seq(idx, k))
        .collect()
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
            k: u32 
        },

        PatternToNum {
            #[structopt(short)]
            pattern: String
        },

        NumToPattern {
            #[structopt(short)]
            num: usize,
            #[structopt(short)]
            k: u32 
        },

        RevComp {
            #[structopt(short)]
            text: String
        },

        PatternMatch {
            #[structopt(short)]
            text: String,
            #[structopt(short="f")]
            is_file: bool,
            #[structopt(short)]
            pattern: String
        },

        Clumps {
            #[structopt(short)]
            text: String,
            #[structopt(short="f")]
            is_file: bool,
            #[structopt(short)]
            k: u32,
            #[structopt(short)]
            l: u32,
            #[structopt(short="h")]
            t: u32
        },
    }
}

use std::path::PathBuf;
use structopt::StructOpt;

fn main() {
    match cli::Command::from_args() {

        cli::Command::CountPattern{text, pattern} => println!("{:?}",
            pattern_count(&seqs::from_str(&text), &seqs::from_str(&pattern))
        ),

        cli::Command::FrequentWords{text, k} => println!("{:?}",
            //frequent_sort(&seqs::from_str(&text), k).iter().map(|seq| seqs::to_str(&seq)).collect::<Vec<_>>()
            frequent_faster(&seqs::from_str(&text), k).iter().map(|seq| seqs::to_str(&seq)).collect::<Vec<_>>()
        ),

        cli::Command::PatternToNum{pattern} => println!("{:?}", 
            seq_to_num(&seqs::from_str(&pattern))
        ),

        cli::Command::NumToPattern{num, k} => println!("{:?}",
            &seqs::to_str(&num_to_seq(num, k))
        ),

        cli::Command::RevComp{text} => println!("{:?}",
            &seqs::to_str(&seqs::revc(&seqs::from_str(&text)))
        ),

        cli::Command::PatternMatch{text, is_file, pattern} => { 
            let seq = match is_file {
                false => seqs::from_str(&text),
                true => seqs::from_file(&PathBuf::from(&text))
            };

            for ind in pattern_match(&seq, &seqs::from_str(&pattern)) {
                print!("{} ", ind.to_string());
            };
        },

        cli::Command::Clumps{text, is_file, k, t, l} => {
            let seq = match is_file {
                false => seqs::from_str(&text),
                true => seqs::from_file(&PathBuf::from(&text))
            };
            println!("{:?}",
                clumps_hash(&seq, k, t, l).iter().map(|seq| seqs::to_str(&seq)).collect::<Vec<_>>()
            )
        },
    }
}
