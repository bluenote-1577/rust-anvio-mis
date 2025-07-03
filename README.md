This is a rust reimplentation of the anvio script for finding misassemblies (https://anvio.org/help/main/programs/anvi-script-find-misassemblies/) by Trigodet et al (https://www.biorxiv.org/content/10.1101/2025.04.22.649783v2). 

Should be quite a bit faster and have near identical results.

### Requirements

- The [rust programming language](https://www.rust-lang.org/)

### Install and usage

```sh
git clone https://github.com/bluenote-1577/rust-anvio-mis
cd rust-anvio-mis
cargo install --path .
bam_error_detector mapping.bam my_file_prefix
```
