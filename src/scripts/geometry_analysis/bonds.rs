use std::fs;
use std::collections::HashMap;
use lazy_static::lazy_static;

const BOND_THRES: f64 = 1.2;

lazy_static! {
    static ref COV_RADS: HashMap<String, f64> = {
        let mut m = HashMap::new();
        m.insert("H".to_string(), 0.37);
        m.insert("C".to_string(), 0.77);
        m.insert("O".to_string(), 0.73);
        m.insert("N".to_string(), 0.75);
        m.insert("F".to_string(), 0.71);
        m.insert("P".to_string(), 1.10);
        m.insert("S".to_string(), 1.03);
        m.insert("Cl".to_string(), 0.99);
        m.insert("Br".to_string(), 1.14);
        m.insert("I".to_string(), 1.33);
        m.insert("He".to_string(), 0.30);
        m.insert("Ne".to_string(), 0.84);
        m.insert("Ar".to_string(), 1.00);
        m.insert("Li".to_string(), 1.02);
        m.insert("Be".to_string(), 0.27);
        m.insert("B".to_string(), 0.88);
        m.insert("Na".to_string(), 1.02);
        m.insert("Mg".to_string(), 0.72);
        m.insert("Al".to_string(), 1.30);
        m.insert("Si".to_string(), 1.18);
        m.insert("K".to_string(), 1.38);
        m.insert("Ca".to_string(), 1.00);
        m.insert("Sc".to_string(), 0.75);
        m.insert("Ti".to_string(), 0.86);
        m.insert("V".to_string(), 0.79);
        m.insert("Cr".to_string(), 0.73);
        m.insert("Mn".to_string(), 0.67);
        m.insert("Fe".to_string(), 0.61);
        m.insert("Co".to_string(), 0.64);
        m.insert("Ni".to_string(), 0.55);
        m.insert("Cu".to_string(), 0.46);
        m.insert("Zn".to_string(), 0.60);
        m.insert("Ga".to_string(), 1.22);
        m.insert("Ge".to_string(), 1.22);
        m.insert("As".to_string(), 1.22);
        m.insert("Se".to_string(), 1.17);
        m.insert("Kr".to_string(), 1.03);
        m.insert("X".to_string(), 0.00);
        m
    };
}


fn get_file_string_array(file_name: &str) -> Vec<Vec<String>> {
    let file_contents = fs::read_to_string(file_name).expect("Could not read file");
    let file_lines: Vec<String> = file_contents.lines().map(|x| x.to_string()).collect();

    let file_string_array = file_lines.iter()
        .map(|x| x.split_whitespace().map(|x| x.to_string()).collect())
        .collect();
    file_string_array
}
fn get_geom(xyz_file_name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let xyz_array = get_file_string_array(xyz_file_name);
    let n_atoms:u32  = xyz_array[0][0].parse().unwrap_or_else(|_| 0);

    // let coords = vec![vec![0.0; 3]; n_atoms as usize];

    // let at_types :Vec<&String> = xyz_array.iter().skip(2).map(|x| &x[0]).collect();
    let at_types :Vec<String> = xyz_array.iter().skip(2)
        .map(|x| x[0].clone()).collect();
    let coords :Vec<Vec<f64>> = xyz_array.iter().skip(2)
        .map(|x| x.iter().skip(1)
            .map(|y| y.parse().unwrap_or_else(|_| 0.0)).collect())
        .collect();
    let geom = (at_types, coords);
    geom
}

fn get_r12(coords1: &Vec<f64>, coords2: &Vec<f64>) -> f64 {
    let mut r2 = 0.0;
    for p in 0..3 {
        r2 += (coords1[p] - coords2[p]).powi(2);
    }
    r2.sqrt()
}

fn get_bond_graph(geom: &(Vec<String>, Vec<Vec<f64>>)) -> Vec<Vec<i64>> {
    let at_types = &geom.0;
    let coords = &geom.1;
    let mut bond_graph: Vec<Vec<i64>> = vec![vec![]; at_types.len()];

    for i in 0..at_types.len(){
        let covrad1 = COV_RADS[&at_types[i]];
        for j in i+1..at_types.len(){
            let covrad2 = COV_RADS[&at_types[j]];
            let thresh = BOND_THRES * (covrad1 + covrad2);
            let r12 = get_r12(&coords[i], &coords[j]);
            if r12 < thresh {
                bond_graph[i].push(j as i64);
                bond_graph[j].push(i as i64);
            }
        }
    }
    bond_graph
}

fn get_bonds(geom: &(Vec<String>, Vec<Vec<f64>>), bond_graph: &Vec<Vec<i64>>) -> Vec<(i64, i64, f64)> {
    let at_types = &geom.0;
    let coords = &geom.1;
    let mut bonds : Vec<(i64, i64, f64)> = vec![];
    for i in 0..at_types.len() {
        for a in 0..bond_graph[i].len() {
            let j = bond_graph[i][a];
            if (i as i64) < j {
                let r12 = get_r12(&coords[i], &coords[j as usize]);
                bonds.push((i as i64, j as i64, r12));
            }
        }
    }
    bonds
}

pub fn bonds(xyz_file_name: &str) {
    let geom = get_geom(xyz_file_name);
    let bond_graph = get_bond_graph(&geom);

    let bonds = get_bonds(&geom, &bond_graph);

    for bond in bonds {
        println!("{} {} {}", bond.0, bond.1, bond.2);
    }
}