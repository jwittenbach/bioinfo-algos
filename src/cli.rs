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