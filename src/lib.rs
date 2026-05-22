use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn nj_from_matrix(input: &Path, output: &mut dyn Write) -> Result<()> {
    let file = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty input".into()))?
        .map_err(RsomicsError::Io)?;
    let names: Vec<&str> = header.split('\t').collect();
    let n = names.len();

    let mut dist = vec![vec![0.0f64; n]; n];
    for (i, line) in lines.enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let vals: Vec<&str> = line.split('\t').collect();
        if vals.len() != n {
            return Err(RsomicsError::InvalidInput(format!(
                "row {i} has {} columns, expected {n}",
                vals.len()
            )));
        }
        for (j, v) in vals.iter().enumerate() {
            dist[i][j] = v
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad float: {e}")))?;
        }
    }

    let newick = neighbor_joining(&names, &dist);
    writeln!(output, "{newick}").map_err(RsomicsError::Io)?;
    Ok(())
}

// Canonical neighbor-joining (full Q-matrix scan). The distance matrix is a flat
// row-major buffer sized for all 2n-1 nodes up front, so the O(N²) inner loops are
// cache-friendly and joins append in place — same tree as the textbook algorithm,
// no per-iteration allocation.
#[allow(clippy::cast_precision_loss, clippy::many_single_char_names)]
fn neighbor_joining(names: &[&str], dist: &[Vec<f64>]) -> String {
    let n = names.len();
    if n <= 1 {
        return names.first().map_or(String::new(), |s| format!("{s};"));
    }

    let max = 2 * n;
    let mut d = vec![0.0f64; max * max];
    for (i, row) in dist.iter().enumerate() {
        let base = i * max;
        for (j, &v) in row.iter().enumerate() {
            d[base + j] = v;
        }
    }

    let mut labels: Vec<String> = Vec::with_capacity(max);
    labels.extend(names.iter().map(|s| (*s).to_string()));
    let mut active = vec![false; max];
    active[..n].fill(true);
    let mut count = n;

    let mut alive: Vec<usize> = Vec::with_capacity(max);
    let mut r = vec![0.0f64; max];

    loop {
        alive.clear();
        alive.extend((0..count).filter(|&i| active[i]));
        let m = alive.len();
        if m <= 2 {
            break;
        }

        let inv = 1.0 / (m - 2) as f64;
        for &i in &alive {
            let base = i * max;
            let mut s = 0.0;
            for &j in &alive {
                s += d[base + j];
            }
            r[i] = s * inv;
        }

        let mut min_val = f64::INFINITY;
        let mut min_i = alive[0];
        let mut min_j = alive[1];
        for (ai, &i) in alive.iter().enumerate() {
            let base = i * max;
            let ri = r[i];
            for &j in &alive[ai + 1..] {
                let q = d[base + j] - ri - r[j];
                if q < min_val {
                    min_val = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let dist_ij = d[min_i * max + min_j];
        let branch_i = 0.5 * dist_ij + 0.5 * (r[min_i] - r[min_j]);
        let branch_j = dist_ij - branch_i;
        labels.push(format!(
            "({}:{branch_i:.6},{}:{branch_j:.6})",
            labels[min_i], labels[min_j]
        ));

        let nu = count;
        let nb = nu * max;
        for &k in &alive {
            if k == min_i || k == min_j {
                continue;
            }
            let dk = 0.5 * (d[k * max + min_i] + d[k * max + min_j] - dist_ij);
            d[nb + k] = dk;
            d[k * max + nu] = dk;
        }
        active[min_i] = false;
        active[min_j] = false;
        active[nu] = true;
        count += 1;
    }

    match alive.as_slice() {
        [i, j] => format!(
            "({}:{:.6},{}:{:.6});",
            labels[*i],
            d[i * max + j] / 2.0,
            labels[*j],
            d[i * max + j] / 2.0
        ),
        [i] => format!("{};", labels[*i]),
        _ => String::from("();"),
    }
}
