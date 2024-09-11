<h1 align="center">abPOA Rust Bindings</h1>
<h2 align="center">Adaptive-band partial order alignment in Rust</h2>

<p>&nbsp;</p>

## Installation

### Cargo package

TODO

### Building from source

#### Rust compiler

The minimum supported Rust version is 1.80.

#### Building abPOA

1. Clone the repository. 

   ```bash
   git clone https://github.com/broadinstitute/abpoa-rs
   ```
2. Move into the directory. 

   ```bash
   cd abpoa-rs
   ```
3. Build using `cargo`. We enable a flag to ensure the compiler uses all features of your machine's CPU. 
   To maximize portability of the binary, however, remove the `RUSTFLAGS="..."` part.

   ```bash
   RUSTFLAGS="-C target-cpu=native" cargo build --release
   ```
   
## Supported features

- [x] Global, local, and semi-global alignment
- [x] Configuring linear, gap-affine, and convex aligment penalties
- [x] Compututing one or more consensus sequences
- [x] Generating row-column MSA FASTA output
- [x] Importing a POA graph from a FASTA file

Features not yet supported:

- [ ] "Strand ambiguous" alignment (on the roadmap)
- [ ] Guide-tree supported alignment
- [ ] Minimizer-based seeding and alignment
   
## Usage

```rust
// Configure the alignment parameters
let aln_params = AlignmentParametersBuilder::new()
    .alignment_mode(AlignmentMode::Global)
    .gap_affine_penalties(0, 4, 6, 2)
    .verbosity(Verbosity::None)
    .build();

// Create a new emoty POA graph
let mut graph = Graph::new(&aln_params);

let test_seq: Vec<&[u8]> = vec![
    b"ACGTGTACAGTTGAC",
    b"AGGTACACGTTAC",
    b"AGTGTCACGTTGAC",
    b"ACGTGTACATTGAC",
];

// Align and add each sequence to the graph
for (i, seq) in test_seq.iter().enumerate() {
    let weights = vec![1; seq.len()];
    graph
        .align_and_add_sequence(&aln_params, seq, &weights, format!("seq{}", i + 1).as_bytes())
        .unwrap();
}

// Compute the row-column MSA output
graph.generate_rc_msa();
let msa = graph.get_msa();
assert_eq!(msa.len(), 4);

let truth = [
    b"ACGTGTACAGTTGAC",
    b"A--GGTACACGTTAC",
    b"A-GTGTCACGTTGAC",
    b"ACGTGTACA-TTGAC",
];

for (i, seq) in msa.sequences().iter().enumerate() {
    // The sequence in the MSA object is coded, so we use `reverse_seq` to convert it to ASCII.
    let ascii = aln_params.reverse_seq(seq);
    assert_eq!(&ascii, truth[i]);
}

// Generate consensus
graph.generate_consensus(ConsensusAlgorithm::HeaviestBundle);

let consensus = graph.get_consensus().unwrap();
let ascii = aln_params.reverse_seq(consensus.sequence());

eprintln!("Consensus: {}", std::str::from_utf8(ascii).unwrap());
```
