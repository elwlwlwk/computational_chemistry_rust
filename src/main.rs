mod scripts;

use scripts::{geometry_analysis::bonds::bonds};
use crate::scripts::geometry_analysis::angles::angles;
use crate::scripts::geometry_analysis::torsions::torsions;

fn main() {
    torsions("src/geom/xyz/ethane.xyz")
}
