use std::error::Error;
use std::ffi::c_void;
use std::fmt::{self, Display};
use std::marker::PhantomData;
use std::path::Path;
use std::sync::LazyLock;

use abpoa_sys::ffi;

pub static LOG_TABLE65536: LazyLock<&'static [u8]> = LazyLock::new(|| unsafe {
    let ptr = ffi::AB_LOG_TABLE_65536.as_ptr();
    std::slice::from_raw_parts(ptr, 65536)
});

pub static BIT_TABLE16: LazyLock<&'static [u8]> = LazyLock::new(|| unsafe {
    let ptr = ffi::AB_BIT_TABLE_16.as_ptr();
    std::slice::from_raw_parts(ptr, 65536)
});

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
        unsafe {
            (*abpoa_params).set_out_msa(1);
            (*abpoa_params).set_out_cons(1);
        }
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

    pub fn gap_affine_penalties(
        self,
        match_score: i32,
        mismatch: i32,
        gap_open: i32,
        gap_ext: i32,
    ) -> Self {
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

    pub fn convex_penalties(
        self,
        match_score: i32,
        mismatch: i32,
        gap_open1: i32,
        gap_ext1: i32,
        gap_open2: i32,
        gap_ext2: i32,
    ) -> Self {
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
        unsafe { (*self.abpoa_params).set_amb_strand(v as u8) }

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
            (*self.abpoa_params).m = matrix_size as i32;
            (*self.abpoa_params).mat = libc::realloc(
                (*self.abpoa_params).mat as *mut libc::c_void,
                matrix_size * matrix_size * std::mem::size_of::<i32>(),
            ) as *mut i32;
        }

        self
    }

    /// Set the filename of an existing POA graph to load.
    ///
    /// The file type can be both GFA or a FASTA with an existing MSA.
    pub fn existing_graph_fname(self, fname: &Path) -> Self {
        let mut fname_cstr = Vec::from(fname.as_os_str().as_encoded_bytes());
        fname_cstr.push(0); // Add null-byte to make it a C-style string

        unsafe {
            (*self.abpoa_params).incr_fn = fname_cstr.as_ptr() as *mut i8;
        }

        // Don't free the C string, as it's now owned by the C struct
        std::mem::forget(fname_cstr);

        self
    }

    /// Use the read IDs in the FASTA or GFA instead of newly set numbers.
    pub fn use_read_ids(self, v: bool) -> Self {
        unsafe {
            (*self.abpoa_params).set_use_read_ids(v as u8);
        }

        self
    }

    /// Construct the final alignment parameters structure.
    ///
    /// Runs some post processing on the paramaters, and initializes additional scoring matrices.
    pub fn build(mut self) -> AlignmentParameters {
        unsafe { ffi::abpoa_post_set_para(self.abpoa_params) }

        let table = unsafe {
            if (*self.abpoa_params).m > 5 {
                std::slice::from_raw_parts(ffi::AB_AA26_TABLE.as_ptr(), 256)
            } else {
                std::slice::from_raw_parts(ffi::AB_NT4_TABLE.as_ptr(), 256)
            }
        };

        let rev_table = unsafe {
            if (*self.abpoa_params).m > 5 {
                std::slice::from_raw_parts(ffi::AB_AA256_TABLE.as_ptr(), 256)
            } else {
                std::slice::from_raw_parts(ffi::AB_NT256_TABLE.as_ptr(), 256)
            }
        };

        LazyLock::force(&LOG_TABLE65536);
        LazyLock::force(&BIT_TABLE16);

        let to_return = AlignmentParameters {
            abpoa_params: self.abpoa_params,
            table,
            rev_table,
        };

        // Set our current parameters to null so we don't double free them (see Drop impl below)
        self.abpoa_params = std::ptr::null_mut();

        to_return
    }

    fn output_msa(self, v: bool) -> Self {
        unsafe {
            (*self.abpoa_params).set_out_msa(v as u8);
        }

        self
    }

    fn output_consensus(self, v: bool) -> Self {
        unsafe {
            (*self.abpoa_params).set_out_cons(v as u8);
        }

        self
    }

    fn consensus_algorithm(self, algorithm: ConsensusAlgorithm) -> Self {
        unsafe {
            (*self.abpoa_params).cons_algrm = algorithm as i32;
        }

        self
    }

    fn max_number_of_consensus_sequences(self, n: i32) -> Self {
        unsafe {
            (*self.abpoa_params).max_n_cons = n;
        }

        self
    }

    fn min_consensus_frequency(self, f: f64) -> Self {
        unsafe {
            (*self.abpoa_params).min_freq = f;
        }

        self
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
    table: &'static [u8],
    rev_table: &'static [u8],
}

impl AlignmentParameters {
    pub fn alphabet_size(&self) -> usize {
        unsafe { (*self.abpoa_params).m as usize }
    }

    pub fn uses_amino_acids(&self) -> bool {
        self.alphabet_size() > 5
    }

    pub fn transform_seq(&self, seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|&x| self.table[x as usize]).collect()
    }

    pub fn reverse_seq(&self, transformed: &[u8]) -> Vec<u8> {
        transformed
            .iter()
            .map(|&x| self.rev_table[x as usize])
            .collect()
    }

    pub fn transform_base(&self, base: u8) -> u8 {
        self.table[base as usize]
    }

    pub fn reverse_base(&self, code: u8) -> u8 {
        self.rev_table[code as usize]
    }
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
        AlignmentParametersBuilder::new().build()
    }
}

/// An enum representing various error conditions that can occur during alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AbpoaError {
    /// The alphabet used in the sequence does not match the graph
    InvalidAlphabet,
    
    /// The alignment input is incorrect
    InvalidInput,
}

impl Error for AbpoaError {}

impl fmt::Display for AbpoaError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AbpoaError::InvalidAlphabet => write!(
                f,
                "The alphabet used in the sequence does not match the graph"
            ),
            AbpoaError::InvalidInput => write!(f, "The alignment input is incorrect"),
        }
    }
}

/// The multiple sequence alignment POA graph
pub struct Graph {
    graph_impl: *mut ffi::abpoa_t,
    use_amino_acids: bool,
}

impl Graph {
    pub fn new(aln_params: &AlignmentParameters) -> Self {
        let graph_impl = unsafe { ffi::abpoa_init() };
        unsafe { ffi::abpoa_reset(graph_impl, aln_params.abpoa_params, 1024) }

        Graph {
            graph_impl,
            use_amino_acids: aln_params.uses_amino_acids(),
        }
    }

    /// Load a graph from a file.
    ///
    /// The function always succeeds. If the file could not be read, it will return an empty graph. Existing
    /// read IDs in the FASTA or GFA will be retained.
    pub fn from_file(fname: &Path, use_amino_acids: bool) -> Self {
        let graph_impl = unsafe { ffi::abpoa_init() };

        let aln_params = AlignmentParametersBuilder::new()
            .existing_graph_fname(fname)
            .use_amino_acids(use_amino_acids)
            .build();

        unsafe {
            ffi::abpoa_restore_graph(graph_impl, aln_params.abpoa_params);
        }

        Graph {
            graph_impl,
            use_amino_acids,
        }
    }

    /// Load a graph from a file.
    ///
    /// The function always succeeds. If the file could not be read, it will return an empty graph.
    /// This function will ignore read IDs in the FASTA or GFA, and use internally generated read IDs. To
    /// retain the existing read IDs, use [`Graph::from_file`].
    pub fn from_file_with_read_ids(fname: &Path, use_amino_acids: bool) -> Self {
        let graph_impl = unsafe { ffi::abpoa_init() };

        let aln_params = AlignmentParametersBuilder::new()
            .existing_graph_fname(fname)
            .use_read_ids(true)
            .use_amino_acids(use_amino_acids)
            .build();

        unsafe {
            ffi::abpoa_restore_graph(graph_impl, aln_params.abpoa_params);
        }

        Graph {
            graph_impl,
            use_amino_acids,
        }
    }

    pub fn sequences(&self) -> Sequences<'_> {
        Sequences::new(self)
    }

    fn get_graph_ptr(&self) -> *const ffi::abpoa_graph_t {
        unsafe { (*self.graph_impl).abg }
    }

    fn get_graph_ptr_mut(&mut self) -> *mut ffi::abpoa_graph_t {
        unsafe { (*self.graph_impl).abg }
    }

    fn get_abs_ptr(&self) -> *const ffi::abpoa_seq_t {
        unsafe { (*self.graph_impl).abs }
    }

    fn get_abs_ptr_mut(&mut self) -> *mut ffi::abpoa_seq_t {
        unsafe { (*self.graph_impl).abs }
    }

    fn get_cons_ptr(&self) -> *const ffi::abpoa_cons_t {
        unsafe { (*self.graph_impl).abc }
    }

    pub fn num_nodes(&self) -> usize {
        unsafe { (*(*self.graph_impl).abg).node_n as usize }
    }

    pub fn num_sequences(&self) -> usize {
        unsafe { (*(*self.graph_impl).abs).n_seq as usize }
    }

    pub fn get_node(&self, node_id: i32) -> NodeRef<'_> {
        unsafe {
            assert!(
                (node_id >= 0) && (node_id as usize) < self.num_nodes(),
                "Node ID out of bounds"
            );

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

    pub fn has_consensus(&self) -> bool {
        unsafe { (*self.get_graph_ptr()).is_called_cons() > 0 && !self.get_cons_ptr().is_null() }
    }
    
    fn align_sequence_coded(
        &mut self,
        aln_params: &AlignmentParameters,
        sequence: &[u8],
    ) -> Result<AlignmentResult, AbpoaError> {
        if aln_params.uses_amino_acids() != self.use_amino_acids {
            return Err(AbpoaError::InvalidAlphabet);
        }

        // Perform alignment
        let mut result = AlignmentResult::new();
        unsafe {
            ffi::abpoa_align_sequence_to_graph(
                self.graph_impl,
                aln_params.abpoa_params,
                sequence.as_ptr() as *mut u8,
                sequence.len() as i32,
                result.as_mut_ptr(),
            );
        }

        Ok(result)
    }
    
    pub fn align_sequence(
        &mut self,
        aln_params: &AlignmentParameters,
        sequence: &[u8],
    ) -> Result<AlignmentResult, AbpoaError> {
        let transformed_seq = aln_params.transform_seq(sequence);
        
        self.align_sequence_coded(aln_params, &transformed_seq)
    }
    
    pub fn add_alignment(
        &mut self,
        aln_params: &AlignmentParameters,
        sequence: &[u8],
        weights: &[i32],
        name: &[u8],
        aln_result: &AlignmentResult,
    ) {
        let (num_existing_seq, num_new) = self.prepare_for_new_sequences(&[name]);

        unsafe {
            ffi::abpoa_add_graph_alignment(
                self.graph_impl,
                aln_params.abpoa_params,
                sequence.as_ptr() as *mut u8,
                weights.as_ptr() as *mut i32,
                sequence.len() as i32,
                std::ptr::null_mut::<i32>(),
                aln_result.result_impl,
                num_existing_seq as i32,
                (num_existing_seq + num_new) as i32,
                1,
            );
        }
    }

    pub fn align_and_add_sequence(
        &mut self,
        aln_params: &AlignmentParameters,
        sequence: &[u8],
        weights: &[i32],
        name: &[u8],
    ) -> Result<AlignmentResult, AbpoaError> {
        let transformed_seq = aln_params.transform_seq(sequence);
        let result = self.align_sequence_coded(aln_params, &transformed_seq)?;
        
        self.add_alignment(aln_params, sequence, weights, name, &result);

        Ok(result)
    }
    
    pub fn align_and_add_multiple<S, W, N>(
        &mut self,
        aln_params: &AlignmentParameters,
        sequences: &[S],
        weights: &[W],
        names: &[N],
    ) -> Result<Vec<AlignmentResult>, AbpoaError>
    where
        S: AsRef<[u8]>,
        W: AsRef<[i32]>,
        N: AsRef<[u8]>,
    {
        if sequences.len() != weights.len() || sequences.len() != names.len() {
            return Err(AbpoaError::InvalidInput);
        }
        
        let transformed_seqs: Vec<_> = sequences.iter()
            .map(|seq| aln_params.transform_seq(seq.as_ref()))
            .collect();
        
        let (num_existing_seq, num_new) = self.prepare_for_new_sequences(names);
        
        let mut all_results = Vec::with_capacity(sequences.len());
        for (i, (seq, w)) in transformed_seqs.iter()
            .zip(weights.iter())
            .enumerate() 
        {
            let result = self.align_sequence_coded(aln_params, seq)?;
            let w = w.as_ref();
            
            unsafe {
                ffi::abpoa_add_graph_alignment(
                    self.graph_impl,
                    aln_params.abpoa_params,
                    seq.as_ptr() as *mut u8,
                    w.as_ptr() as *mut i32,
                    seq.len() as i32,
                    std::ptr::null_mut::<i32>(),
                    result.result_impl,
                    (num_existing_seq + i) as i32,
                    (num_existing_seq + num_new) as i32,
                    1,
                );
            }
            
            all_results.push(result);
        }
        
        Ok(all_results)
    }

    fn prepare_for_new_sequences<N: AsRef<[u8]>>(&mut self, names: &[N]) -> (usize, usize) {
        let num_new_sequences = names.len();
        let num_existing_seq = unsafe { (*self.get_abs_ptr()).n_seq as usize };
        
        unsafe {
            (*self.get_abs_ptr_mut()).n_seq += num_new_sequences as i32;
            ffi::abpoa_realloc_seq(self.get_abs_ptr_mut());
        }
        
        // Set new sequence names
        for (i, name) in names.iter().enumerate() {
            unsafe {
                let n = name.as_ref();
                let target = (*self.get_abs_ptr()).name.add(num_existing_seq + i);
                ffi::abpoa_cpy_str(target, n.as_ptr() as *mut i8, n.len() as i32)
            }
        }
        
        (num_existing_seq, num_new_sequences)
    }

    pub fn generate_consensus(&mut self, alg: ConsensusAlgorithm) {
        let cons_aln_params = AlignmentParametersBuilder::new()
            .output_consensus(true)
            .consensus_algorithm(alg)
            .use_amino_acids(self.use_amino_acids)
            .build();

        if alg == ConsensusAlgorithm::Majority {
            self.alloc_node_msa_rank();
        }

        unsafe {
            ffi::abpoa_generate_consensus(self.graph_impl, cons_aln_params.abpoa_params);
        }
    }

    pub fn force_generate_consensus(&mut self, alg: ConsensusAlgorithm) {
        unsafe {
            (*self.get_graph_ptr_mut()).set_is_called_cons(0);
        }

        self.generate_consensus(alg)
    }

    pub fn generate_consensus_multiple(
        &mut self,
        alg: ConsensusAlgorithm,
        max_consensus: usize,
        min_freq: Option<f64>,
    ) {
        let cons_aln_params = AlignmentParametersBuilder::new()
            .output_consensus(true)
            .consensus_algorithm(alg)
            .use_amino_acids(self.use_amino_acids)
            .max_number_of_consensus_sequences(max_consensus as i32)
            .min_consensus_frequency(min_freq.unwrap_or(0.25))
            .build();

        self.alloc_node_msa_rank();
        unsafe {
            ffi::abpoa_generate_consensus(self.graph_impl, cons_aln_params.abpoa_params);
        }
    }

    pub fn force_generate_consensus_multiple(
        &mut self,
        alg: ConsensusAlgorithm,
        max_consensus: usize,
        min_freq: Option<f64>,
    ) {
        unsafe {
            (*self.get_graph_ptr_mut()).set_is_called_cons(0);
        }

        self.generate_consensus_multiple(alg, max_consensus, min_freq)
    }

    pub fn generate_rc_msa(&mut self) {
        let aln_params = AlignmentParametersBuilder::new()
            .output_msa(true)
            .output_consensus(self.has_consensus())
            .use_amino_acids(self.use_amino_acids)
            .build();

        self.alloc_node_msa_rank();

        unsafe {
            ffi::abpoa_generate_rc_msa(self.graph_impl, aln_params.abpoa_params);
        }
    }

    pub fn get_consensus(&self) -> Option<ConsensusData<'_>> {
        if self.has_consensus() {
            Some(ConsensusData::new(self))
        } else {
            None
        }
    }

    pub fn get_msa(&self) -> MSAData<'_> {
        MSAData::new(self)
    }

    fn alloc_node_msa_rank(&mut self) {
        let graph_ptr = self.get_graph_ptr_mut();
        unsafe {
            (*graph_ptr).node_id_to_msa_rank = libc::realloc(
                (*graph_ptr).node_id_to_msa_rank as *mut libc::c_void,
                (*graph_ptr).index_rank_m as usize * std::mem::size_of::<i32>(),
            ) as *mut i32;

            if (*graph_ptr).node_id_to_msa_rank.is_null() {
                panic!("Failed to allocate memory for node_id_to_msa_rank");
            }
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

/// Reference to a node in the graph.
///
/// This struct is initialized by calling the `get_node` method on a `Graph` object.
pub struct NodeRef<'a> {
    node_data: *const ffi::abpoa_node_t,
    graph: PhantomData<&'a Graph>,
}

impl<'a> NodeRef<'a> {
    fn new(_: &'a Graph, node_data: *const ffi::abpoa_node_t) -> Self {
        NodeRef {
            node_data,
            graph: PhantomData,
        }
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

/// View of the sequences in the POA graph
pub struct Sequences<'graph> {
    seq: Vec<AbpoaStr>,
    names: Vec<AbpoaStr>,
    comments: Vec<AbpoaStr>,
    qualities: Vec<AbpoaStr>,
    is_rc: Vec<bool>,
    num_sequences: usize,
    graph: PhantomData<&'graph Graph>,
}

impl<'graph> Sequences<'graph> {
    fn new(graph: &'graph Graph) -> Self {
        let seq_ptr = graph.get_abs_ptr();
        let num_sequences = unsafe { (*seq_ptr).n_seq as usize };
        unsafe {
            Sequences {
                seq: if (*seq_ptr).seq.is_null() {
                    Vec::default()
                } else {
                    std::slice::from_raw_parts((*seq_ptr).seq, num_sequences)
                        .iter()
                        .map(AbpoaStr::from_ref)
                        .collect()
                },
                names: if (*seq_ptr).name.is_null() {
                    Vec::default()
                } else {
                    std::slice::from_raw_parts((*seq_ptr).name, num_sequences)
                        .iter()
                        .map(AbpoaStr::from_ref)
                        .collect()
                },
                comments: if (*seq_ptr).comment.is_null() {
                    Vec::default()
                } else {
                    std::slice::from_raw_parts((*seq_ptr).comment, num_sequences)
                        .iter()
                        .map(AbpoaStr::from_ref)
                        .collect()
                },
                qualities: if (*seq_ptr).qual.is_null() {
                    Vec::default()
                } else {
                    std::slice::from_raw_parts((*seq_ptr).qual, num_sequences)
                        .iter()
                        .map(AbpoaStr::from_ref)
                        .collect()
                },
                is_rc: if (*seq_ptr).is_rc.is_null() {
                    Vec::default()
                } else {
                    std::slice::from_raw_parts((*seq_ptr).is_rc, num_sequences)
                        .iter()
                        .map(|v| *v != 0)
                        .collect()
                },
                num_sequences: graph.num_sequences(),
                graph: PhantomData,
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.num_sequences == 0
    }

    pub fn len(&self) -> usize {
        self.num_sequences
    }

    pub fn sequences(&self) -> &[AbpoaStr] {
        &self.seq
    }

    pub fn names(&self) -> &[AbpoaStr] {
        &self.names
    }

    pub fn comments(&self) -> &[AbpoaStr] {
        &self.comments
    }

    pub fn qualities(&self) -> &[AbpoaStr] {
        &self.qualities
    }

    pub fn iter(&self) -> SequenceIter {
        SequenceIter::new(self)
    }
}

/// Iterator over the sequences in the POA graph
pub struct SequenceIter<'a> {
    sequences: &'a Sequences<'a>,
    index: usize,
}

impl<'a> SequenceIter<'a> {
    fn new(sequences: &'a Sequences) -> Self {
        SequenceIter {
            sequences,
            index: 0,
        }
    }
}

impl<'a> Iterator for SequenceIter<'a> {
    type Item = Sequence<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.sequences.len() {
            let seq = Sequence {
                sequences: self.sequences,
                index: self.index,
            };
            self.index += 1;
            Some(seq)
        } else {
            None
        }
    }
}

impl<'a> ExactSizeIterator for SequenceIter<'a> {
    fn len(&self) -> usize {
        self.sequences.len()
    }
}

/// A reference to a particular sequence in the graph
pub struct Sequence<'a> {
    sequences: &'a Sequences<'a>,
    index: usize,
}

impl<'a> Sequence<'a> {
    pub fn sequence(&self) -> &[u8] {
        self.sequences
            .seq
            .get(self.index)
            .map(|v| v.as_ref())
            .unwrap_or(b"")
    }

    pub fn name(&self) -> &[u8] {
        self.sequences
            .names
            .get(self.index)
            .map(|v| v.as_ref())
            .unwrap_or(b"")
    }

    pub fn comment(&self) -> &[u8] {
        self.sequences
            .comments
            .get(self.index)
            .map(|v| v.as_ref())
            .unwrap_or(b"")
    }

    pub fn quality(&self) -> &[u8] {
        self.sequences
            .qualities
            .get(self.index)
            .map(|v| v.as_ref())
            .unwrap_or(b"")
    }

    pub fn is_rc(&self) -> bool {
        self.sequences
            .is_rc
            .get(self.index)
            .copied()
            .unwrap_or(false)
    }
}

/// A non-owning reference to an abPOA string.
pub struct AbpoaStr {
    str_impl: *const ffi::abpoa_str_t,
}

impl AbpoaStr {
    fn from_ref(string: &ffi::abpoa_str_t) -> Self {
        AbpoaStr {
            str_impl: std::ptr::addr_of!(*string),
        }
    }

    pub fn is_empty(&self) -> bool {
        if self.str_impl.is_null() {
            return true;
        }

        self.len() == 0 || unsafe { (*self.str_impl).s.is_null() }
    }

    pub fn len(&self) -> usize {
        unsafe { (*self.str_impl).l as usize }
    }

    pub fn capacity(&self) -> usize {
        unsafe { (*self.str_impl).m as usize }
    }
}

impl std::ops::Index<usize> for AbpoaStr {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.len() {
            panic!("index out of bounds");
        }

        unsafe {
            let ptr = (*self.str_impl).s as *const u8;
            ptr.add(index).as_ref().unwrap()
        }
    }
}

impl PartialEq for AbpoaStr {
    fn eq(&self, other: &Self) -> bool {
        *self.as_ref() == *other.as_ref()
    }
}

impl Eq for AbpoaStr {}

impl PartialEq<str> for AbpoaStr {
    fn eq(&self, other: &str) -> bool {
        *self.as_ref() == *other.as_bytes()
    }
}

impl PartialEq<&str> for AbpoaStr {
    fn eq(&self, other: &&str) -> bool {
        *self.as_ref() == *other.as_bytes()
    }
}

impl PartialEq<[u8]> for AbpoaStr {
    fn eq(&self, other: &[u8]) -> bool {
        *self.as_ref() == *other
    }
}

impl PartialEq<&[u8]> for AbpoaStr {
    fn eq(&self, other: &&[u8]) -> bool {
        *self.as_ref() == **other
    }
}

impl<T, const N: usize> PartialEq<[T; N]> for AbpoaStr
where
    T: PartialEq<u8>,
{
    fn eq(&self, other: &[T; N]) -> bool {
        *other == self.as_ref()
    }
}

impl<T, const N: usize> PartialEq<&[T; N]> for AbpoaStr
where
    T: PartialEq<u8>,
{
    fn eq(&self, other: &&[T; N]) -> bool {
        **other == self.as_ref()
    }
}

impl Display for AbpoaStr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_empty() {
            return write!(f, "");
        }

        let string = unsafe {
            std::str::from_utf8(std::slice::from_raw_parts(
                (*self.str_impl).s as *const u8,
                self.len(),
            ))
            .unwrap_or("[error: invalid utf-8]")
        };
        write!(f, "{}", string)
    }
}

impl fmt::Debug for AbpoaStr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "AbpoaStr({})", self)
    }
}

impl AsRef<[u8]> for AbpoaStr {
    fn as_ref(&self) -> &[u8] {
        if self.is_empty() {
            b""
        } else {
            unsafe { std::slice::from_raw_parts((*self.str_impl).s as *const u8, self.len()) }
        }
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

    pub fn print_alignment(
        &self,
        aln_params: &AlignmentParameters,
        graph: &Graph,
        sequence: &[u8],
    ) {
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
                    let is_match = aln_params.reverse_base(graph.get_node(node_id).get_base())
                        == sequence[query_id];
                    graph_bases.push(aln_params.reverse_base(graph.get_node(node_id).get_base()));
                    aln_symbols.push(if is_match { b'|' } else { b'.' });
                    qry_bases.push(sequence[query_id]);
                    eprintln!("{len}{op}");
                }
                CigarOp::Deletion => {
                    graph_bases.push(graph.get_node(node_id).get_base());
                    aln_symbols.push(b' ');
                    qry_bases.push(b'-');
                    eprintln!("{len}D");
                }
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
            unsafe { libc::free(self.result_impl.graph_cigar as *mut c_void) }
        }
    }
}

pub struct ConsensusData<'a> {
    consensus_data_impl: *const ffi::abpoa_cons_t,
    seq: Vec<&'a [u8]>,
}

impl<'a> ConsensusData<'a> {
    fn new(graph: &'a Graph) -> Self {
        let cons_ptr = graph.get_cons_ptr();
        let num_seq = unsafe { (*cons_ptr).n_cons as usize };

        let seqs = (0..num_seq)
            .map(|i| {
                let len = unsafe { (*cons_ptr).cons_len.add(i) };
                unsafe { std::slice::from_raw_parts(*(*cons_ptr).cons_base.add(i), *len as usize) }
            })
            .collect();

        ConsensusData {
            consensus_data_impl: graph.get_cons_ptr(),
            seq: seqs,
        }
    }

    pub fn len(&self) -> usize {
        unsafe { (*self.consensus_data_impl).n_cons as usize }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn consensus_lengths(&self) -> &[i32] {
        unsafe { std::slice::from_raw_parts((*self.consensus_data_impl).cons_len, self.len()) }
    }

    pub fn sequences(&self) -> &[&[u8]] {
        &self.seq
    }
}

pub struct MSAData<'a> {
    consensus_data_impl: *const ffi::abpoa_cons_t,
    alignments: Vec<&'a [u8]>,
}

impl<'a> MSAData<'a> {
    fn new(graph: &'a Graph) -> Self {
        let cons_ptr = graph.get_cons_ptr();
        let num_seq = unsafe { (*cons_ptr).n_seq as usize + (*cons_ptr).n_cons as usize };
        let msa_len = unsafe { (*cons_ptr).msa_len as usize };

        let alignments = (0..num_seq)
            .map(|i| unsafe { std::slice::from_raw_parts(*(*cons_ptr).msa_base.add(i), msa_len) })
            .collect();

        MSAData {
            consensus_data_impl: graph.get_cons_ptr(),
            alignments,
        }
    }

    pub fn includes_consesus(&self) -> bool {
        unsafe { (*self.consensus_data_impl).n_cons > 0 }
    }

    pub fn len(&self) -> usize {
        unsafe { (*self.consensus_data_impl).n_seq as usize }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn sequences(&self) -> &[&[u8]] {
        &self.alignments
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

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
            .verbosity(Verbosity::None)
            .build();

        let mut graph = Graph::new(&aln_params);
        assert_eq!(graph.num_nodes(), 2); // Includes the special source and sink nodes

        let sequence = b"CGTACGTACTGACGTACGATCGTACTGACGTCGTCA";
        let weights = vec![1; sequence.len()];

        let result = graph
            .align_and_add_sequence(&aln_params, sequence, &weights, b"seq1")
            .unwrap();
        assert_eq!(result.get_best_score(), 0);

        for (i, orig_base) in sequence.iter().enumerate() {
            let node = graph.get_node(i as i32 + 2);
            let base = aln_params.reverse_base(node.get_base());
            assert_eq!(base, *orig_base);
        }
    }

    #[test]
    fn test_graph_from_file() {
        let test_fname = PathBuf::from("tests/data/small_test.truth.fa");
        let graph = Graph::from_file(&test_fname, false);

        assert_eq!(graph.num_sequences(), 3);
        assert_eq!(graph.sequences().names(), &[b"seq1", b"seq2", b"seq3"]);
    }

    #[test]
    fn test_alignment_to_graph() {
        let aln_params = AlignmentParametersBuilder::new()
            .alignment_mode(AlignmentMode::Global)
            .gap_affine_penalties(0, 4, 6, 2)
            .verbosity(Verbosity::Debug)
            .build();

        let mut graph = Graph::new(&aln_params);
        assert_eq!(graph.num_nodes(), 2); // Includes the special source and sink nodes
        assert_eq!(graph.num_sequences(), 0);

        let sequence = b"CGTACGTACTGACGTACGATCGTACTGACGTCGTCA";
        let weights = vec![1; sequence.len()];

        let result = graph
            .align_and_add_sequence(&aln_params, sequence, &weights, b"seq1")
            .unwrap();
        assert_eq!(result.get_best_score(), 0);
        assert_eq!(graph.num_nodes(), 2 + sequence.len());
        assert_eq!(graph.num_sequences(), 1);

        let sequence2 = b"CGTACGTACTGACGTTTGATCGTACTGACGTCGTCA";
        let weights2 = vec![1; sequence2.len()];

        let result2 = graph
            .align_and_add_sequence(&aln_params, sequence2, &weights2, b"seq2")
            .unwrap();

        assert_eq!(result2.get_best_score(), -8);
        assert_eq!(result2.get_num_matches(), sequence2.len() - 2);
        assert_eq!(graph.num_nodes(), 2 + sequence.len() + 2);
    }
    
    #[test]
    fn test_align_multiple() {
        let aln_params = AlignmentParametersBuilder::new()
            .alignment_mode(AlignmentMode::Global)
            .gap_affine_penalties(0, 4, 6, 2)
            .verbosity(Verbosity::None)
            .build();

        let mut graph = Graph::new(&aln_params);

        let test_seq: Vec<&[u8]> = vec![
            b"ACGTGTACAGTTGAC",
            b"AGGTACACGTTAC",
            b"AGTGTCACGTTGAC",
            b"ACGTGTACATTGAC",
        ];
        
        let weights: Vec<_> = test_seq.iter()
            .map(|seq| vec![1i32; seq.len()])
            .collect();
        
        let names: Vec<_> = (1..=test_seq.len())
            .map(|i| format!("seq{}", i))
            .collect();
        
        let _ = graph.align_and_add_multiple(&aln_params, &test_seq, weights.as_slice(), &names)
            .unwrap();
        
        assert_eq!(graph.num_sequences(), 4);
    }

    #[test]
    fn test_consensus_generation() {
        let aln_params = AlignmentParametersBuilder::new()
            .alignment_mode(AlignmentMode::Global)
            .gap_affine_penalties(0, 4, 6, 2)
            .verbosity(Verbosity::None)
            .build();

        let mut graph =
            Graph::from_file(&PathBuf::from("tests/data/test_from_abpoa.truth.fa"), false);

        graph.generate_consensus(ConsensusAlgorithm::HeaviestBundle);

        let consensus = graph.get_consensus().unwrap();
        let truth = aln_params.transform_seq(b"ACGTGTACAGTTGAC");

        assert_eq!(consensus.len(), 1);
        assert_eq!(consensus.sequences(), &[&truth]);
    }

    #[test]
    fn test_msa_generation() {
        let aln_params = AlignmentParametersBuilder::new()
            .alignment_mode(AlignmentMode::Global)
            .gap_affine_penalties(0, 4, 6, 2)
            .verbosity(Verbosity::None)
            .build();

        let mut graph = Graph::new(&aln_params);

        let test_seq: Vec<&[u8]> = vec![
            b"ACGTGTACAGTTGAC",
            b"AGGTACACGTTAC",
            b"AGTGTCACGTTGAC",
            b"ACGTGTACATTGAC",
        ];

        for (i, seq) in test_seq.iter().enumerate() {
            let weights = vec![1; seq.len()];
            graph
                .align_and_add_sequence(&aln_params, seq, &weights, format!("{}", i + 1).as_bytes())
                .unwrap();
        }

        eprintln!("Generate MSA");
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
            let ascii = aln_params.reverse_seq(seq);
            assert_eq!(&ascii, truth[i]);
        }
    }
}
