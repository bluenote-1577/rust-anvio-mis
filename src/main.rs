// BAM long read assembly error detector
// Converted from Python to Rust
// This script identifies potential errors in long read assemblies by analyzing:
// 1. Regions with clipping at both ends
// 2. Regions with zero coverage

use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::Write;
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::CigarStringView;
use clap::{Arg, Command};

#[derive(Debug)]
struct ContigData {
    coverage: Vec<u32>,
    clipping: HashMap<usize, u32>,
    length: usize,
}

impl ContigData {
    fn new(length: usize) -> Self {
        let coverage = vec![0;length];
        
        ContigData {
            coverage,
            clipping: HashMap::new(),
            length,
        }
    }
    
    fn add_clipping(&mut self, pos: usize) {
        *self.clipping.entry(pos).or_insert(0) += 1;
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command line arguments
    let matches = Command::new("bam_error_detector")
        .version("1.0")
        .author("Converted from anvi'o script")
        .about("Identifies potential errors in long read assemblies using mapping information")
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Input BAM file of long reads mapped to an assembly made from these reads"))
        .arg(Arg::new("output_prefix")
            .required(true)
            .help("Prefix for output files"))
        .arg(Arg::new("min_dist_to_end")
            .long("min-dist-to-end")
            .default_value("100")
            .help("Minimum distance from contig ends to consider"))
        .arg(Arg::new("clipping_ratio")
            .long("clipping-ratio")
            .default_value("1.0")
            .help("Minimum ratio of clipped reads to total coverage to report"))
        .get_matches();

    let bam_file = matches.get_one::<String>("bam_file").unwrap();
    let output_prefix = matches.get_one::<String>("output_prefix").unwrap();
    let min_dist_to_end: usize = matches.get_one::<String>("min_dist_to_end").unwrap().parse()?;
    let min_clipping_ratio: f64 = matches.get_one::<String>("clipping_ratio").unwrap().parse()?;
    let just_do_it = true;

    if !just_do_it {
        eprintln!("This script ONLY makes sense if you are using a BAM file that was made from");
        eprintln!("mapping long read onto an assembly made with the SAME long reads.");
        eprintln!("If you are positive that you did JUST that, then re-run this program with");
        eprintln!("--just-do-it flag.");
        return Err("Missing --just-do-it flag".into());
    }

    println!("BAM file: {}", bam_file);
    println!("Length of contig's end to ignore: {}", min_dist_to_end);

    // Open BAM file
    let mut bam = bam::Reader::from_path(bam_file)?;
    
    // Get header information to extract reference names and lengths
    let header = bam.header().clone();
    
    // Create a map of contig names to lengths
    let mut contigs_size = HashMap::new();
    for (i, contig) in header.target_names().iter().enumerate() {
        let contig_name = String::from_utf8_lossy(contig).to_string();
        let contig_len = header.target_len(i as u32).unwrap() as usize;
        contigs_size.insert(contig_name, contig_len);
    }
    
    // The main dictionary to store coverage information
    let mut cov_dict: HashMap<String, ContigData> = HashMap::new();
    
    
    // Read counter
    let mut read_count = 0;
    
    // Process each read in the BAM file
    for result in bam.records() {
        let read = result?;
        read_count += 1;
        
        if read.is_unmapped() {
            continue;
        }
        
        if read_count % 500 == 0 {
            eprint!("\rProcessed {} reads", read_count);
        }
        
        let contig = String::from_utf8_lossy(header.tid2name(read.tid() as u32)).to_string();
        let contig_end = contigs_size.get(&contig).unwrap();
        
        // Initialize the contig data if it doesn't exist
        if !cov_dict.contains_key(&contig) {
            cov_dict.insert(contig.clone(), ContigData::new(*contig_end));
        }
        
        let mut current_pos = read.pos() as usize;
        let cigar = read.cigar();
        
        // Count the number of CIGAR operations
        let mut num_tup = 0;
        let contig_struct = cov_dict.get_mut(&contig).unwrap();
        
        for op in cigar.iter() {
            num_tup += 1;
            
            match op.char() {
                // If mapping (M, =, X), compute coverage, increase current position
                'M' | '=' | 'X' => {
                    for pos in current_pos..(current_pos + op.len() as usize) {
                        contig_struct.coverage[pos] += 1
                    }
                    current_pos += op.len() as usize;
                },
                // If deletion (D), increase current position
                'D' => {
                    current_pos += op.len() as usize;
                },
                // If clipping (S, H), then +1 clipping at position
                'S' | 'H' => {
                    if num_tup == 1 {
                        // If at start of contig, skip
                        if current_pos != 0 {
                            contig_struct.add_clipping(current_pos);
                        }
                    } else if current_pos != *contig_end {
                        contig_struct.add_clipping(current_pos.saturating_sub(1));
                    }
                },
                _ => {}
            }
        }
    }
    println!("\nRead processing complete");

    // Write clipping output file
    let clipping_output = format!("{}-clipping.txt", output_prefix);
    println!("Output file: {}", clipping_output);
    
    let mut clipping_file = File::create(&clipping_output)?;
    writeln!(clipping_file, "contig\tlength\tpos\trelative_pos\tcov\tclipping\tclipping_ratio")?;
    
    for (contig, data) in &cov_dict {
        let contig_length = data.length;
        
        for (&pos, &clipping) in &data.clipping {
            let cov = data.coverage[pos];
            let clipping_ratio = clipping as f64 / cov as f64;
            let relative_pos = pos as f64 / contig_length as f64;
            
            if clipping_ratio >= min_clipping_ratio && 
               pos > min_dist_to_end && 
               contig_length - pos > min_dist_to_end {
                writeln!(
                    clipping_file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    contig, contig_length, pos, relative_pos, cov, clipping, clipping_ratio
                )?;
            }
        }
    }
    
    // Write zero coverage output file
    let zero_output = format!("{}-zero_cov.txt", output_prefix);
    println!("Output file: {}", zero_output);
    
    let mut zero_file = File::create(&zero_output)?;
    writeln!(zero_file, "contig\tlength\trange\trange_size")?;
    
    for (contig, data) in &cov_dict {
        let contig_length = data.length;
        let mut in_window = false;
        let mut window_start = 0;
        
        for pos in 0..contig_length {
            if data.coverage[pos] == 0 && !in_window {
                window_start = pos;
                in_window = true;
                write!(zero_file, "{}\t{}\t{}-", contig, contig_length, window_start)?;
            } else if data.coverage[pos] > 0 && in_window {
                let window_end = pos;
                let window_length = window_end - window_start;
                in_window = false;
                writeln!(zero_file, "{}\t{}", window_end, window_length)?;
            }
            
            // If end of contig
            if data.coverage[0] == 0 && pos == contig_length - 1 {
                if in_window {
                    let window_end = pos + 1;
                    let window_length = window_end - window_start;
                    writeln!(zero_file, "{}\t{}", window_end, window_length)?;
                } else {
                    let window_start = pos;
                    let window_end = pos + 1;
                    let window_length = window_end - window_start;
                    writeln!(
                        zero_file,
                        "{}\t{}\t{}-{}\t{}",
                        contig, contig_length, window_start, window_end, window_length
                    )?;
                }
            }
        }

        if in_window{
            writeln!(zero_file, "");
        }
    }
    
    println!("Analysis complete!");
    Ok(())
}
