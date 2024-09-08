mod scripts;

use scripts::{geometry_analysis::bonds::bonds};

fn main() {
    bonds("src/geom/xyz/ethane.xyz")
}
