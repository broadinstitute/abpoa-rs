use std::ffi::c_void;

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

impl Default for Graph {
    fn default() -> Self {
        let graph_impl = unsafe { ffi::abpoa_init() };
        Graph { graph_impl }
    }
}

impl Graph {
    pub fn new() -> Self {
        Self::default()
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
    
    pub fn get_num_aligned(&self) -> i32 {
        self.result_impl.n_aln_bases
    }
    
    pub fn get_num_matches(&self) -> i32 {
        self.result_impl.n_matched_bases
    }
}

impl Drop for AlignmentResult {
    fn drop(&mut self) {
        if self.result_impl.n_cigar > 0 {
            unsafe { libc::free(self.result_impl.graph_cigar as *mut c_void)}
        }
    }
}

/// Align a single sequence to a abPOA graph
pub fn align_to_graph(graph: &mut Graph, aln_params: &AlignmentParameters, sequence: &[u8]) -> AlignmentResult {
    let mut result = AlignmentResult::new();
    unsafe {
        ffi::abpoa_align_sequence_to_graph(graph.graph_impl, aln_params.abpoa_params, sequence.as_ptr() as *mut u8, sequence.len() as i32, result.as_mut_ptr());
    }
    result
}


/// Perform a multiple sequence alignment
pub fn msa() {
    
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_align_params_builder() {
        let params = AlignmentParametersBuilder::new()
            .ambiguous_strand(true)
            .band_width(1)
            .band_f(1.0)
            .use_amino_acids(true)
            .consensus_algorithm(ConsensusAlgorithm::HeaviestBundle)
            .max_number_of_consensus_sequences(1)
            .min_consensus_frequency(0.5)
            .build();
        
        assert_eq!(params.abpoa_params.is_null(), false);
        ass
    }
}