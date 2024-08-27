mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

#[repr(u32)]
pub enum AlignmentMode {
    Global = ffi::ABPOA_GLOBAL_MODE,
    Local = ffi::ABPOA_LOCAL_MODE,
    Extend = ffi::ABPOA_EXTEND_MODE,
}

#[repr(u32)]
pub enum AlignmentType {
    Linear = ffi::ABPOA_LINEAR_GAP,
    Affine = ffi::ABPOA_AFFINE_GAP,
    Convex = ffi::ABPOA_CONVEX_GAP,
}

pub struct AlignmentParameters {
    abpoa_params: *mut ffi::abpoa_para_t,
}

impl AlignmentParameters {
    pub fn new() -> Self {
        let abpoa_params = unsafe { ffi::abpoa_init_para() };
        AlignmentParameters { abpoa_params }
    }
}

impl Drop for AlignmentParameters {
    fn drop(&mut self) {
        unsafe { ffi::abpoa_free_para(self.abpoa_params) };
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
        unsafe { ffi::abpoa_free(self.graph_impl) };
    }
}
