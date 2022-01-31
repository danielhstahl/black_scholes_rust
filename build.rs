extern crate cc;

fn main() {
    cc::Build::new()
        .cpp(true)
        .file("./letsberational/erf_cody.cpp")
        .file("./letsberational/rationalcubic.cpp")
        .file("./letsberational/normaldistribution.cpp")
        .file("./letsberational/lets_be_rational.cpp")
        .compile("letsberational.a");
}
