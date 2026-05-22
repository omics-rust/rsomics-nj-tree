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

// Neighbor-joining via rapidNJ-style bound-based search (Simonsen, Mailund &
// Pedersen 2008). Per row the distances are kept sorted ascending; the minimum
// Q-value Q(i,j)=d(i,j)-u_i-u_j is found by scanning each row until
// d(i,j)-u_i-u_max ≥ current best (no later entry can beat it). This yields the
// identical tree to the canonical full O(N²)-per-step scan while skipping most of
// the matrix. Distances never change for surviving pairs, so the sorted rows stay
// valid; inactive columns are skipped lazily and each new node gets its own sorted
// row (which carries every pair that involves it).
#[allow(
    clippy::cast_precision_loss,
    clippy::many_single_char_names,
    clippy::needless_range_loop
)]
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

    // Raw row sums over the active set (d(i,i)=0 contributes nothing).
    let mut r = vec![0.0f64; max];
    for i in 0..n {
        let base = i * max;
        r[i] = (0..n).map(|j| d[base + j]).sum();
    }

    // Per-row distances sorted ascending (column j != i).
    let mut srow: Vec<Vec<(f64, u32)>> = vec![Vec::new(); max];
    for i in 0..n {
        let base = i * max;
        let mut v: Vec<(f64, u32)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (d[base + j], j as u32))
            .collect();
        v.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        srow[i] = v;
    }

    let mut count = n;
    let mut m = n;

    while m > 2 {
        let inv = 1.0 / (m - 2) as f64;
        let mut u_max = f64::NEG_INFINITY;
        for i in 0..count {
            if active[i] {
                let u = r[i] * inv;
                if u > u_max {
                    u_max = u;
                }
            }
        }

        let mut q_min = f64::INFINITY;
        let mut mi = usize::MAX;
        let mut mj = usize::MAX;
        for i in 0..count {
            if !active[i] {
                continue;
            }
            let u_i = r[i] * inv;
            for &(dij, jr) in &srow[i] {
                if dij - u_i - u_max >= q_min {
                    break;
                }
                let j = jr as usize;
                if !active[j] {
                    continue;
                }
                let q = dij - u_i - r[j] * inv;
                if q < q_min {
                    q_min = q;
                    mi = i;
                    mj = j;
                }
            }
        }

        let dist_ij = d[mi * max + mj];
        let branch_i = 0.5 * dist_ij + 0.5 * (r[mi] - r[mj]) * inv;
        let branch_j = dist_ij - branch_i;
        labels.push(format!(
            "({}:{branch_i:.6},{}:{branch_j:.6})",
            labels[mi], labels[mj]
        ));

        let nu = count;
        let nb = nu * max;
        let mut r_nu = 0.0;
        for k in 0..count {
            if !active[k] || k == mi || k == mj {
                continue;
            }
            let dk = 0.5 * (d[k * max + mi] + d[k * max + mj] - dist_ij);
            r[k] = r[k] - d[k * max + mi] - d[k * max + mj] + dk;
            d[nb + k] = dk;
            d[k * max + nu] = dk;
            r_nu += dk;
        }
        r[nu] = r_nu;

        active[mi] = false;
        active[mj] = false;
        active[nu] = true;

        let mut v: Vec<(f64, u32)> = (0..count)
            .filter(|&k| active[k] && k != nu)
            .map(|k| (d[nb + k], k as u32))
            .collect();
        v.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        srow[nu] = v;

        count += 1;
        m -= 1;
    }

    let rest: Vec<usize> = (0..count).filter(|&i| active[i]).collect();
    match rest.as_slice() {
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
