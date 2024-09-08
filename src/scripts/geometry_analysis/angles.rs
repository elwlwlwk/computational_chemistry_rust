use crate::scripts::geometry_analysis::bonds::{get_bond_graph, get_bonds, get_geom, get_r12, print_bonds, print_geom};

fn get_u12(coords1: &Vec<f64>, coords2: &&Vec<f64>) -> Vec<f64> {
    let r12 = get_r12(&coords1, &coords2);
    let u12 = coords1.iter().zip(coords2.iter()).map(|(x, y)| (x - y) / r12).collect();
    u12
}

fn get_udp(uvec1: &Vec<f64>, uvec2: &Vec<f64>) -> f64 {
    let udp = uvec1.iter().zip(uvec2.iter()).map(|(x, y)| x * y).sum();
    udp
}

fn get_a123(xyz1: &Vec<f64>, xyz2: &Vec<f64>, xyz3: &Vec<f64>) -> f64 {
    let u12 = get_u12(xyz1, &xyz2);
    let u32 = get_u12(xyz3, &xyz2);
    let cos12 = get_udp(&u12, &u32);
    let a123 = cos12.acos();
    a123
}

pub fn get_angles(geom: &(Vec<String>, Vec<Vec<f64>>), bond_graph: &Vec<Vec<i64>>) -> Vec<(i64, i64, i64, f64)> {
    let at_types = &geom.0;
    let coords = &geom.1;

    let mut angles : Vec<(i64, i64, i64, f64)> = vec![];
    for j in 0..at_types.len() {
        let n_jbonds = bond_graph[j].len();
        for a in 0..n_jbonds {
            let i = bond_graph[j][a];
            for b in a+1..n_jbonds {
                let k = bond_graph[j][b];
                let a123: f64 = get_a123(&coords[i as usize], &coords[j], &coords[k as usize]);
                angles.append(&mut vec![(i, j as i64, k, a123)]);
            }
        }
    }

    angles
}

fn print_angles(geom: &(Vec<String>, Vec<Vec<f64>>), angles: &Vec<(i64, i64, i64, f64)>) {
    let at_types = &geom.0;

    println!("{} angle(s) found (degrees)", angles.len());
    for (n1, n2, n3, a123) in angles {
        let nstr = format!("{}-{}-{}", n1 + 1, n2 + 1, n3 + 1);
        let tstr = format!("({}-{}-{})", at_types[*n1 as usize], at_types[*n2 as usize], at_types[*n3 as usize]);
        println!("{: <6} {: <10} {: <10}", nstr, tstr, a123.to_degrees());
    }
    println!();
}

pub fn angles(xyz_file_name: &str) {
    let geom = get_geom(xyz_file_name);
    let bond_graph = get_bond_graph(&geom);

    let bonds = get_bonds(&geom, &bond_graph);
    let angles = get_angles(&geom, &bond_graph);

    print_geom(&geom, &String::from("inital geometry"));
    print_bonds(&geom, &bonds);
    print_angles(&geom, &angles);
}