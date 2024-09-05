use std::{ffi::c_void, fmt, marker::PhantomData};

mod ffi {
    #![allow(non_upper_case_globals)]
    #![allow(non_camel_case_types)]
    #![allow(non_snake_case)]
    #![allow(dead_code)]
    
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

#[derive(PartialEq, Eq, Debug, Copy, Clone)]
#[repr(u32)]
pub enum AlignmentMode {
    Global = ffi::ABPOA_GLOBAL_MODE,
    Local = ffi::ABPOA_LOCAL_MODE,
    Extend = ffi::ABPOA_EXTEND_MODE,
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
#[repr(u32)]
pub enum ConsensusAlgorithm {
    HeaviestBundle = ffi::ABPOA_HB,
    Majority = ffi::ABPOA_MC,
}


#[derive(PartialEq, Eq, Debug, Clone, Copy, PartialOrd, Ord)]
#[repr(u32)]
pub enum Verbosity {
    None = ffi::ABPOA_NONE_VERBOSE,
    Info = ffi::ABPOA_INFO_VERBOSE,
    Debug = ffi::ABPOA_DEBUG_VERBOSE,
    LongDebug = ffi::ABPOA_LONG_DEBUG_VERBOSE,
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
#[repr(u32)]
pub enum CigarOp {
    Match = ffi::ABPOA_CMATCH,
    Insertion = ffi::ABPOA_CINS,
    Deletion = ffi::ABPOA_CDEL,
    Mismatch = ffi::ABPOA_CDIFF,
    SoftClip = ffi::ABPOA_CSOFT_CLIP,
    HardClip = ffi::ABPOA_CHARD_CLIP,
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CigarOp::Match => write!(f, "M"),
            CigarOp::Insertion => write!(f, "I"),
            CigarOp::Deletion => write!(f, "D"),
            CigarOp::Mismatch => write!(f, "X"),
            CigarOp::SoftClip => write!(f, "S"),
            CigarOp::HardClip => write!(f, "H"),
        }
    }
}

impl From<u64> for CigarOp {
    fn from(mut value: u64) -> Self {
        value &= 0xf;
        
        match value {
            0 => CigarOp::Match,
            1 => CigarOp::Insertion,
            2 => CigarOp::Deletion,
            3 => CigarOp::Mismatch,
            4 => CigarOp::SoftClip,
            5 => CigarOp::HardClip,
            _ => panic!("Invalid CIGAR operation"),
        }
    }
}

pub struct AlignmentParametersBuilder {
    abpoa_params: *mut ffi::abpoa_para_t,
}

impl AlignmentParametersBuilder {
    pub fn new() -> Self {
        let abpoa_params = unsafe { ffi::abpoa_init_para() };
        AlignmentParametersBuilder { abpoa_params }
    }
    
    pub fn alignment_mode(self, mode: AlignmentMode) -> Self {
        unsafe {
            (*self.abpoa_params).align_mode = mode as i32;
        }
        
        self
    }
    
    pub fn verbosity(self, mode: Verbosity) -> Self {
        unsafe {
            (*self.abpoa_params).verbose = mode as i32;
        }
        
        self
    }
    
    pub fn gap_linear_penalties(self, match_score: i32, mismatch: i32, gap: i32) -> Self {
        unsafe {
            (*self.abpoa_params).match_ = match_score;
            (*self.abpoa_params).mismatch = mismatch;
            (*self.abpoa_params).gap_open1 = 0;
            (*self.abpoa_params).gap_ext1 = gap;
            (*self.abpoa_params).gap_open2 = 0;
            (*self.abpoa_params).gap_ext2 = 0;
        }
        
        self
    }
    
    pub fn gap_affine_penalties(self, match_score: i32, mismatch: i32, gap_open: i32, gap_ext: i32) -> Self {
        unsafe {
            (*self.abpoa_params).match_ = match_score;
            (*self.abpoa_params).mismatch = mismatch;
            (*self.abpoa_params).gap_open1 = gap_open;
            (*self.abpoa_params).gap_ext1 = gap_ext;
            (*self.abpoa_params).gap_open2 = 0;
            (*self.abpoa_params).gap_ext2 = 0;
        }
        
        self
    }
    
    pub fn convex_penalties(self, match_score: i32, mismatch: i32, gap_open1: i32, gap_ext1: i32, gap_open2: i32, gap_ext2: i32) -> Self {
        unsafe {
            (*self.abpoa_params).match_ = match_score;
            (*self.abpoa_params).mismatch = mismatch;
            (*self.abpoa_params).gap_open1 = gap_open1;
            (*self.abpoa_params).gap_ext1 = gap_ext1;
            (*self.abpoa_params).gap_open2 = gap_open2;
            (*self.abpoa_params).gap_ext2 = gap_ext2;
        }
        
        self
    }
    
    pub fn ambiguous_strand(self, v: bool) -> Self {
        unsafe {
            (*self.abpoa_params).set_amb_strand(v as u8)
        }
        
        self
    }
    
    // Adaptive band width. Set to 0 to disable.
    pub fn band_width(self, w: i32) -> Self {
        unsafe {
            (*self.abpoa_params).wb = w;
        }
        
        self
    }
    
    pub fn band_f(self, f: f32) -> Self {
        unsafe {
            (*self.abpoa_params).wf = f;
        }
        
        self
    }
    
    pub fn use_amino_acids(self, v: bool) -> Self {
        let matrix_size: usize = if v { 27 } else { 5 };
        unsafe {
            libc::free((*self.abpoa_params).mat as *mut c_void);
            
            (*self.abpoa_params).m = matrix_size as i32;
            (*self.abpoa_params).mat = libc::malloc(matrix_size * matrix_size * std::mem::size_of::<i32>()) as *mut i32;
        }
        
        self
    }
    
    pub fn consensus_algorithm(self, algorithm: ConsensusAlgorithm) -> Self {
        unsafe {
            (*self.abpoa_params).cons_algrm = algorithm as i32;
        }
        
        self
    }
    
    pub fn max_number_of_consensus_sequences(self, n: i32) -> Self {
        unsafe {
            (*self.abpoa_params).max_n_cons = n;
        }
        
        self
    }
    
    pub fn min_consensus_frequency(self, f: f64) -> Self {
        unsafe {
            (*self.abpoa_params).min_freq = f;
        }
        
        self
    }
    
    /// Construct the final alignment parameters structure.
    /// 
    /// Runs some post processing on the paramaters, and initializes additional scoring matrices.
    pub fn build(mut self) -> AlignmentParameters {
        unsafe {
            ffi::abpoa_post_set_para(self.abpoa_params)
        }
        
        let to_return = AlignmentParameters { abpoa_params: self.abpoa_params };
        
        // Set our current parameters to null so we don't double free them (see Drop impl below)
        self.abpoa_params = std::ptr::null_mut();
        
        to_return
    }
}

impl Default for AlignmentParametersBuilder {
    fn default() -> Self {
        AlignmentParametersBuilder::new()
    }
}

impl Drop for AlignmentParametersBuilder {
    fn drop(&mut self) {
        if self.abpoa_params.is_null() {
            return;
        }
        
        unsafe { ffi::abpoa_free_para(self.abpoa_params) };
    }
}

pub struct AlignmentParameters {
    abpoa_params: *mut ffi::abpoa_para_t,
}
    
impl Drop for AlignmentParameters {
    fn drop(&mut self) {
        if self.abpoa_params.is_null() {
            return;
        }
        
        unsafe { ffi::abpoa_free_para(self.abpoa_params) };
    }
}

impl Default for AlignmentParameters {
    fn default() -> Self {
        AlignmentParametersBuilder::new()
            .build()
    }
}


/// The multiple sequence alignment POA graph
pub struct Graph {
    graph_impl: *mut ffi::abpoa_t,
}

impl Graph {
    pub fn new(aln_params: &AlignmentParameters) -> Self {
        let graph_impl = unsafe { ffi::abpoa_init() };
        unsafe { ffi::abpoa_reset(graph_impl, aln_params.abpoa_params, 1024)}
        
        Graph { graph_impl }
    }
    
    fn get_graph_ptr(&self) -> *const ffi::abpoa_graph_t {
        unsafe { (*self.graph_impl).abg } 
    }
    
    fn get_abs_ptr(&mut self) -> *mut ffi::abpoa_seq_t {
        unsafe { (*self.graph_impl).abs }
    }
    
    pub fn num_nodes(&self) -> usize {
        unsafe { (*(*self.graph_impl).abg).node_n as usize }
    }
    
    pub fn num_sequences(&self) -> usize {
        unsafe { (*(*self.graph_impl).abs).n_seq as usize }
    }
    
    pub fn get_node(&self, node_id: i32) -> NodeRef<'_> {
        unsafe {
            assert!((node_id >= 0) && (node_id as usize) < self.num_nodes(), "Node ID out of bounds");
            
            let node_data = (*self.get_graph_ptr()).node.add(node_id as usize);
            NodeRef::new(self, node_data)
        }
    }
    
    pub fn get_node_by_index(&self, node_index: usize) -> NodeRef<'_> {
        assert!(node_index < self.num_nodes(), "Node index out of bounds");
        unsafe {
            let node_id = *(*self.get_graph_ptr()).index_to_node_id.add(node_index);
            self.get_node(node_id)
        }
    }
    
    pub fn align_and_add_sequence(&mut self, aln_params: &AlignmentParameters, sequence: &[u8], weights: &[i32], name: &[u8]) -> AlignmentResult {
        let num_existing_seq = self.num_sequences();
        
        // Make space for the new sequence
        self.prepare_for_new_sequences(1);
        
        // Set new sequence name
        unsafe {
            let target = (*self.get_abs_ptr()).name.add(num_existing_seq);
            ffi::abpoa_cpy_str(target, name.as_ptr() as *mut i8, name.len() as i32)
        }
        
        // Perform alignment
        let mut result = AlignmentResult::new();
        unsafe {
            ffi::abpoa_align_sequence_to_graph(
                self.graph_impl, 
                aln_params.abpoa_params, 
                sequence.as_ptr() as *mut u8, 
                sequence.len() as i32, 
                result.as_mut_ptr()
            );
        }
        
        unsafe {
            ffi::abpoa_add_graph_alignment(
                self.graph_impl, 
                aln_params.abpoa_params, 
                sequence.as_ptr() as *mut u8, 
                weights.as_ptr() as *mut i32, 
                sequence.len() as i32, 
                std::ptr::null_mut::<i32>(), 
                result.result_impl, 
                num_existing_seq as i32, 
                num_existing_seq as i32 + 1, 
                1
            );
        }
        
        result
    }
    
    fn prepare_for_new_sequences(&mut self, num_new_sequences: usize) {
        unsafe {
            (*self.get_abs_ptr()).n_seq += num_new_sequences as i32;
            ffi::abpoa_realloc_seq(self.get_abs_ptr());
        }
    }
}

impl Drop for Graph {
    fn drop(&mut self) {
        if self.graph_impl.is_null() {
            return;
        }
        
        unsafe { ffi::abpoa_free(self.graph_impl) };
    }
}

pub struct NodeRef<'a> {
    node_data: *const ffi::abpoa_node_t,
    graph: PhantomData<&'a Graph>
}

impl<'a> NodeRef<'a> {
    fn new(_: &'a Graph, node_data: *const ffi::abpoa_node_t) -> Self {
        NodeRef { node_data, graph: PhantomData }
    }
    
    pub fn get_id(&self) -> i32 {
        unsafe { (*self.node_data).node_id }
    }
    
    pub fn get_base(&self) -> u8 {
        unsafe { (*self.node_data).base }
    }
    
    pub fn get_base_char(&self) -> char {
        char::from(self.get_base())
    }
    
    pub fn in_degree(&self) -> usize {
        unsafe { (*self.node_data).in_edge_n as usize }
    }
    
    pub fn out_degree(&self) -> usize {
        unsafe { (*self.node_data).out_edge_n as usize }
    }
    
    pub fn num_aligned_nodes(&self) -> usize {
        unsafe { (*self.node_data).aligned_node_n as usize }
    }
}

/// A sequence-to-graph alignment result
pub struct AlignmentResult {
    result_impl: ffi::abpoa_res_t,
}

impl AlignmentResult {
    pub fn new() -> Self {
        let result_impl = unsafe { std::mem::zeroed() };
        AlignmentResult { result_impl }
    }
    
    pub fn as_mut_ptr(&mut self) -> *mut ffi::abpoa_res_t {
        &mut self.result_impl
    }
    
    pub fn get_best_score(&self) -> i32 {
        self.result_impl.best_score
    }
    
    pub fn get_num_aligned(&self) -> usize {
        self.result_impl.n_aln_bases as usize
    }
    
    pub fn get_num_matches(&self) -> usize {
        self.result_impl.n_matched_bases as usize
    }
    
    pub fn print_alignment(&self, graph: &Graph, sequence: &[u8]) {
        let mut graph_bases = Vec::with_capacity(self.get_num_aligned());
        let mut aln_symbols = Vec::with_capacity(self.get_num_aligned());
        let mut qry_bases = Vec::with_capacity(self.get_num_aligned());
        
        for i in 0..self.result_impl.n_cigar {
            let cigar = unsafe { *self.result_impl.graph_cigar.add(i as usize) };
            let op = CigarOp::from(cigar);
            let node_id = (cigar >> 34) as i32;
            let len = (cigar >> 4) & 0x3fffffff;
            let mut query_id = len as usize;
            
            match op {
                CigarOp::Match | CigarOp::Mismatch => {
                    graph_bases.push(graph.get_node(node_id).get_base());
                    aln_symbols.push(if op == CigarOp::Match { b'|' } else { b'.' });
                    qry_bases.push(sequence[query_id]);
                    eprintln!("{len}{op}");
                },
                CigarOp::Deletion => {
                    graph_bases.push(graph.get_node(node_id).get_base());
                    aln_symbols.push(b' ');
                    qry_bases.push(b'-');
                    eprintln!("{len}D");
                },
                CigarOp::Insertion | CigarOp::SoftClip | CigarOp::HardClip => {
                    query_id = node_id as usize;
                    let start = query_id + 1 - len as usize;
                    eprintln!("Query ID: {query_id}");
                    graph_bases.extend(b"-".repeat(len as usize));
                    aln_symbols.extend(b" ".repeat(len as usize));
                    qry_bases.extend(sequence[start..=query_id].iter());
                    eprintln!("{len}{op}");
                }
            }
        }
        
        eprintln!();
        eprintln!("{}", String::from_utf8_lossy(&graph_bases));
        eprintln!("{}", String::from_utf8_lossy(&aln_symbols));
        eprintln!("{}", String::from_utf8_lossy(&qry_bases));
    }
}

impl Default for AlignmentResult {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for AlignmentResult {
    fn drop(&mut self) {
        if self.result_impl.n_cigar > 0 {
            unsafe { libc::free(self.result_impl.graph_cigar as *mut c_void)}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_align_params_builder() {
        let params = AlignmentParametersBuilder::new()
            .ambiguous_strand(true)
            .band_width(10)
            .band_f(1.0)
            .use_amino_acids(true)
            .consensus_algorithm(ConsensusAlgorithm::HeaviestBundle)
            .max_number_of_consensus_sequences(1)
            .min_consensus_frequency(0.5)
            .build();
        
        assert!(!params.abpoa_params.is_null());
        unsafe {
            assert!((*params.abpoa_params).m > 5);
        }
    }
    
    #[test]
    fn test_node_ref() {
        let aln_params = AlignmentParametersBuilder::new()
            .verbosity(Verbosity::Debug)
            .build();
        
        let mut graph = Graph::new(&aln_params);
        assert_eq!(graph.num_nodes(), 2);  // Includes the special source and sink nodes
        
        let sequence = b"CGTACGTACTGACGTACGATCGTACTGACGTCGTCA";
        let weights = vec![1; sequence.len()];
        
        let result = graph.align_and_add_sequence(&aln_params, sequence, &weights, b"seq1");
        assert_eq!(result.get_best_score(), 0);
        
        for (i, orig_base) in sequence.iter().enumerate() {
            let node = graph.get_node(i as i32 + 2);
            let base = node.get_base();
            assert_eq!(base, *orig_base);
        }
    }
    
    #[test]
    fn test_alignment_to_graph() {
        let aln_params = AlignmentParametersBuilder::new()
            .verbosity(Verbosity::Debug)
            .build();
        
        let mut graph = Graph::new(&aln_params);
        assert_eq!(graph.num_nodes(), 2);  // Includes the special source and sink nodes
        assert_eq!(graph.num_sequences(), 0);
        
        eprintln!("Seq 1");
        let sequence = b"CGTACGTACTGACGTACGATCGTACTGACGTCGTCA";
        let weights = vec![1; sequence.len()];
        
        let result = graph.align_and_add_sequence(&aln_params, sequence, &weights, b"seq1");
        assert_eq!(result.get_best_score(), 0);
        assert_eq!(graph.num_nodes(), 2 + sequence.len());
        assert_eq!(graph.num_sequences(), 1);
        
        eprintln!("Seq 2");
        let sequence2 = b"CGTACGTACTGACGTTTGATCGTACTGACGTCGTCA";
        let weights2 = vec![1; sequence2.len()];
        
        let result2 = graph.align_and_add_sequence(&aln_params, sequence2, &weights2, b"seq2");
        result2.print_alignment(&graph, sequence);
        
        eprintln!("{:?} {:?} {:?}", result2.get_best_score(), result2.get_num_aligned(), result2.get_num_matches());
        assert_eq!(result2.get_best_score(), 114);
        assert_eq!(result2.get_num_matches(), sequence2.len() - 2);
        assert_eq!(graph.num_nodes(), 2 + sequence.len() + 2);
    }
}