use crate::scripts::geometry_analysis::angles::{get_angles, get_u12, get_udp};
use crate::scripts::geometry_analysis::bonds::{get_bond_graph, get_bonds, get_geom, get_r12};

const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
const DEG2RAD: f64 = 1.0 / RAD2DEG;

pub fn get_ucp(uvec1: &Vec<f64>, uvec2: &Vec<f64>) -> Vec<f64> {
    let cos_12 = get_udp(uvec1, uvec2);
    let sin_12 = (1.0 - cos_12.powi(2)).sqrt();
    let ucp = vec![
        (uvec1[1] * uvec2[2] - uvec1[2] * uvec2[1]) / sin_12,
        (uvec1[2] * uvec2[0] - uvec1[0] * uvec2[2]) / sin_12,
        ( uvec1[0] * uvec2[1] - uvec1[1] * uvec2[0]) / sin_12
    ];
    ucp
}

fn get_t1234(coords1: &Vec<f64>, coords2: &Vec<f64>, coords3: &Vec<f64>, coords4: &Vec<f64>) -> f64 {
    let u21 = get_u12(&coords2, &coords1);
    let u23 = get_u12(&coords2, &coords3);
    let u32 = get_u12(&coords3, &coords2);
    let u34 = get_u12(&coords3, &coords4);
    let u21c23 = get_ucp(&u21, &u23);
    let u32c34 = get_ucp(&u32, &u34);
    let dp = get_udp(&u21c23, &u32c34);
    let sign = get_udp(&u21c23, &u34).signum() * -1.0;
    let t1234 = RAD2DEG * sign * dp.acos();
    t1234
}

pub fn get_torsions(geom: &(Vec<String>, Vec<Vec<f64>>), bond_graph: &Vec<Vec<i64>>) -> Vec<(i64, i64, i64, i64, f64)> {
    let at_types = &geom.0;
    let coords = &geom.1;
    let n_atoms = at_types.len();

    let mut torsions : Vec<(i64, i64, i64, i64, f64)> = vec![];
    for j in 0..n_atoms {
        let n_jbonds = bond_graph[j].len();
        for a in 0..n_jbonds {
            let k: usize = bond_graph[j][a] as usize;
            if (k < j) {
                continue;
            }
            let n_kbonds = bond_graph[k].len();
            for b in 0..n_kbonds {
                let i = bond_graph[j][b] as usize;
                if (i == k) {
                    continue;
                }
                for c in 0..n_kbonds {
                    let l = bond_graph[k][c] as usize;
                    if (l == j || l == i) {
                        continue;
                    }
                    let t1234 = get_t1234(&coords[i as usize], &coords[j], &coords[k as usize], &coords[l as usize]);
                    torsions.append(&mut vec![(i as i64, j as i64, k as i64, l as i64, t1234)]);
                }
            }
        }
    }
    torsions
}

fn print_torsions(geom: &(Vec<String>, Vec<Vec<f64>>), torsions: &Vec<(i64, i64, i64, i64, f64)>) {
    let at_types = &geom.0;

    println!("{} torsion(s) found (degrees)", torsions.len());
    for (n1, n2, n3, n4, t1234) in torsions {
        let nstr = format!("{}-{}-{}-{}", n1 + 1, n2 + 1, n3 + 1, n4 + 1);
        let tstr = format!("({}-{}-{}-{})", at_types[*n1 as usize], at_types[*n2 as usize], at_types[*n3 as usize], at_types[*n4 as usize]);
        println!("{: <6} {: <10} {: <10}", nstr, tstr, t1234);
    }
    println!();
}

pub fn torsions(xyz_file_name: &str) {
    let geom = get_geom(xyz_file_name);
    let bond_graph = get_bond_graph(&geom);

    let torsions = get_torsions(&geom, &bond_graph);

    print_torsions(&geom, &torsions);
}

// 9 torsion(s) found (degrees)
// 3-1-2-6          (H-C-C-H)        60.001
// 3-1-2-7          (H-C-C-H)       -59.999
// 3-1-2-8          (H-C-C-H)      -180.000
// 4-1-2-6          (H-C-C-H)       -60.000
// 4-1-2-7          (H-C-C-H)      -180.000
// 4-1-2-8          (H-C-C-H)        59.999
// 5-1-2-6          (H-C-C-H)       180.000
// 5-1-2-7          (H-C-C-H)        60.000
// 5-1-2-8          (H-C-C-H)       -60.000