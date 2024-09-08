mod scripts;

use scripts::{geometry_analysis::bonds::bonds};
use crate::scripts::geometry_analysis::angles::angles;

fn main() {
    angles("src/geom/xyz/ethane.xyz")
}
