use crate::config::CONFIG;

const GRID_COLS: usize = CONFIG.physics.grid_cols;
const GRID_ROWS: usize = CONFIG.physics.grid_rows;

fn first_last_mut<T>(slice: &mut [T]) -> (&mut T, &mut T) {
    let (a, b) = slice.split_first_mut().unwrap();
    (a, b.last_mut().unwrap())
}
fn double_index_mut<T>(slice: &mut [T], i: usize) -> (&mut T, &mut T) {
    assert!(slice.len() > i + 1);
    let (a, b) = slice.split_at_mut(i + 1);
    (&mut a[i], &mut b[0])
}
fn triple_index_mut<T>(slice: &mut [T], i: usize) -> (&mut T, &mut T, &mut T) {
    assert!(slice.len() > i + 2);
    let (a, rest) = slice.split_at_mut(i + 1);
    let (b, c) = rest.split_first_mut().unwrap();
    (&mut a[i], b, &mut c[0])
}

pub fn stream_row_range(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    old_lattice_rows: &[[[f32; 9]; GRID_COLS]],
    new_lattice_rows: &mut [[[f32; 9]; GRID_COLS]],
    scratch_rows: &mut [[[f32; 9]; GRID_COLS]],
    base_row: usize,
    leaf_row_count: usize,
) {
    assert_eq!(old_lattice_rows.len(), new_lattice_rows.len());
    if old_lattice_rows.len() <= leaf_row_count {
        let (row_before, row_after) = first_last_mut(scratch_rows);
        zero_out_new_rows(row_before, row_after, new_lattice_rows);
        if old_lattice_rows.len() == 1 {
            stream_row(
                blocked,
                old_lattice_rows.first().unwrap(),
                [row_before, &mut new_lattice_rows[0], row_after],
                base_row,
            );
        } else {
            let (first, second) = double_index_mut(new_lattice_rows, 0);
            stream_row(
                blocked,
                &old_lattice_rows[0],
                [row_before, first, second],
                base_row,
            );
            for i in 0..(old_lattice_rows.len() - 2) {
                let (a, b, c) = triple_index_mut(new_lattice_rows, i);
                stream_row(
                    blocked,
                    &old_lattice_rows[i + 1],
                    [a, b, c],
                    base_row + i + 1,
                );
            }
            let (secondlast, last) = double_index_mut(new_lattice_rows, new_lattice_rows.len() - 2);
            stream_row(
                blocked,
                old_lattice_rows.last().unwrap(),
                [secondlast, last, row_after],
                base_row + old_lattice_rows.len() - 1,
            );
        }
    } else {
        let mid = old_lattice_rows.len() / 2;
        let (old_left, old_right) = old_lattice_rows.split_at(mid);
        let (new_left, new_right) = new_lattice_rows.split_at_mut(mid);
        let (scratch_left, scratch_right) = scratch_rows.split_at_mut(scratch_rows.len() / 2);
        rayon::join(
            || {
                stream_row_range(
                    blocked,
                    old_left,
                    new_left,
                    scratch_left,
                    base_row,
                    leaf_row_count,
                )
            },
            || {
                stream_row_range(
                    blocked,
                    old_right,
                    new_right,
                    scratch_right,
                    base_row + mid,
                    leaf_row_count,
                )
            },
        );
        stitch_together_using_scratch_rows(
            new_left.last_mut().unwrap(),
            new_right.first_mut().unwrap(),
            scratch_left.last().unwrap(),
            scratch_right.first().unwrap(),
        );
    }
}
fn zero_out_new_rows(
    row_before: &mut [[f32; 9]; GRID_COLS],
    row_after: &mut [[f32; 9]; GRID_COLS],
    new_lattice_rows: &mut [[[f32; 9]; GRID_COLS]],
) {
    row_before.fill([0.0; 9]);
    row_after.fill([0.0; 9]);
    new_lattice_rows.fill([[0.0; 9]; GRID_COLS]);
}
fn stitch_together_using_scratch_rows(
    new_left: &mut [[f32; 9]; GRID_COLS],
    new_right: &mut [[f32; 9]; GRID_COLS],
    scratch_left: &[[f32; 9]; GRID_COLS],
    scratch_right: &[[f32; 9]; GRID_COLS],
) {
    type Col = [[f32; 9]; GRID_COLS];
    type Flat = [f32; 9 * GRID_COLS];
    // SAFETY: Col and Flat have the same representation
    unsafe {
        let new_left: &mut Flat = &mut *(new_left as *mut Col as *mut Flat);
        let new_right: &mut Flat = &mut *(new_right as *mut Col as *mut Flat);
        let scratch_left: &Flat = &*(scratch_left as *const Col as *const Flat);
        let scratch_right: &Flat = &*(scratch_right as *const Col as *const Flat);
        // Only after flattening is this loop autovectorized.
        for i in 0..(9 * GRID_COLS) {
            new_left[i] += scratch_right[i];
            new_right[i] += scratch_left[i];
        }
    }
}
fn stream_row(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    old_row: &[[f32; 9]; GRID_COLS],
    new_rows: [&mut [[f32; 9]; GRID_COLS]; 3],
    r0: usize,
) {
    for (c0, old_cell) in old_row.iter().enumerate() {
        if blocked[r0][c0] {
            continue;
        }
        for old_vel in 0..9 {
            let ((r, c), new_vel) = step_cell_coordinates(blocked, (r0, c0), old_vel);
            new_rows[r.wrapping_sub(r0).wrapping_add(1)][c][new_vel] += old_cell[old_vel];
        }
    }
}

fn step_cell_coordinates(
    blocked: &[[bool; GRID_COLS]; GRID_ROWS],
    (r, c): (usize, usize),
    old_vel: usize,
) -> ((usize, usize), usize) {
    let dr = (old_vel / 3) as isize - 1;
    let dc = (old_vel % 3) as isize - 1;

    // Bounce at upper and lower wall
    let (nr, dr) = {
        let nr = r as isize + dr;
        if nr == -1 {
            (0, 1)
        } else if nr == GRID_ROWS as isize {
            (GRID_ROWS - 1, -1)
        } else {
            (nr as usize, dr)
        }
    };
    // Periodic at left and right wall
    let nc = (c.wrapping_add(dc as usize).wrapping_add(GRID_COLS)) % GRID_COLS;

    // Bounce against blocked cells
    let blocked_both = blocked[nr][nc];
    let blocked_r = blocked[nr][c];
    let blocked_c = blocked[r][nc];
    let (r1, c1, dr1, dc1) =
        if (blocked_r && blocked_c) || (blocked_both && !blocked_r && !blocked_c) {
            (r, c, -dr, -dc)
        } else if !blocked_both {
            (nr, nc, dr, dc)
        } else if blocked_c {
            (nr, c, dr, -dc)
        } else {
            (r, nc, -dr, dc)
        };
    ((r1, c1), (dr1 * 3 + dc1 + 4) as usize)
}
