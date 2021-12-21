use delaunator::{triangulate, Point};
use itertools::Itertools;
use std::fs::File;
use std::io::prelude::*;
use std::{cmp, ops};

#[derive(Clone, Debug, Copy)]
struct Vert {
    x: f64,
    y: f64,
}

#[derive(Debug, Clone, Copy)]
struct Tri(usize, usize, usize);

impl Tri {
    fn new(verts: &Vec<Vert>, indices: [usize; 3]) -> Tri {
        for idx in indices {
            if !(idx < verts.len()) {
                panic!("Cannot build Tri from indices outside of referenced Array");
            }
        }
        Tri {
            0: indices[0],
            1: indices[1],
            2: indices[2],
        }
    }
}

impl cmp::PartialEq<Vert> for Vert {
    fn eq(&self, other: &Vert) -> bool {
        self.x == other.x && self.y == other.y
    }

    fn ne(&self, other: &Vert) -> bool {
        !self.eq(other)
    }
}

impl cmp::PartialOrd<Vert> for Vert {
    fn gt(&self, other: &Vert) -> bool {
        self.y > other.y || (self.y == other.y && self.x > other.x)
    }

    fn lt(&self, other: &Vert) -> bool {
        self.y < other.y
    }

    fn ge(&self, other: &Vert) -> bool {
        self.gt(other) && self.eq(other)
    }

    fn le(&self, other: &Vert) -> bool {
        self.lt(other) && self.eq(other)
    }

    fn partial_cmp(&self, other: &Vert) -> Option<cmp::Ordering> {
        Some(if self < other {
            cmp::Ordering::Less
        } else if self > other {
            cmp::Ordering::Greater
        } else {
            cmp::Ordering::Equal
        })
    }
}

impl ops::Add<Vert> for Vert {
    type Output = Vert;
    fn add(self, _rhs: Vert) -> Vert {
        Vert {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
        }
    }
}

impl ops::Sub<Vert> for Vert {
    type Output = Vert;
    fn sub(self, _rhs: Vert) -> Vert {
        Vert {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
        }
    }
}

impl Vert {
    fn verts_2_points(verts: Vec<Vert>) -> Vec<Point> {
        let mut res = Vec::<Point>::new();
        res.reserve(verts.len());
        for vert in verts {
            res.push(Point {
                x: vert.x,
                y: vert.y,
            });
        }
        res
    }

    fn scale_ret(self, fac: f64) -> Vert {
        Vert {
            x: self.x * fac,
            y: self.y * fac,
        }
    }

    fn scale(&mut self, fac: f64) {
        self.x *= fac;
        self.y *= fac;
    }

    fn normalize_ret(&mut self) -> Vert {
        let mag: f64 = (self.x.powi(2) + self.y.powi(2)).sqrt();
        self.scale_ret(1.0 / mag)
    }

    fn normalize(&mut self) {
        let mag: f64 = (self.x.powi(2) + self.y.powi(2)).sqrt();
        self.scale(1.0 / mag);
    }

    fn as_del_point(&self) -> Point {
        Point {
            x: self.x,
            y: self.y,
        }
    }
}

#[derive(Debug, Clone)]
struct SimCell {
    tri: Tri,
    neighbors: [usize; 3],
    isboundary: bool,
    center: Vert,
    value: f64,
}

impl SimCell {
    fn vert_in(&self, p: &Vert, verts: &Vec<Vert>) -> bool {
        let e0 = verts[self.tri.0] - verts[self.tri.1];
        let e1 = verts[self.tri.0] - verts[self.tri.2];
        let c0 = verts[self.tri.0];

        let t1: f64 = (e0.x * (p.y - c0.y) - e0.y * (p.x - c0.x)) / (e1.y * e0.x - e1.x * e0.y);
        let t0: f64 = (p.x - c0.x - t1 * e1.x) / e0.x;

        t0 + t1 < 1.0 && t0 > 0.0 && t1 > 0.0
    }
}

#[derive(Debug, Clone)]
struct Domain {
    verts: Vec<Vert>,
    border: Vec<Vert>,
    cells: Vec<SimCell>,
    area: Vec<Tri>,
}

impl Domain {
    #[allow(dead_code, unused_variables)]
    fn is_inside(border: &Vec<Vert>, p: &Vert) -> bool {
        let filtered = border.windows(2).filter(|line_seg| -> bool {
                let c0 = line_seg[0];
                let c1 = line_seg[1];
                if c0.y < c1.y {
                    p.y >= c0.y && p.y < c1.y
                } else if c0.y > c1.y {
                    p.y >= c1.y && p.y < c0.y
                } else {
                    false
                }
        }); 
        println!("{:?}: {:?}", p, filtered.size_hint().1.or_else(|| Some(1)));
        for line_seg in border.windows(2).filter(|line_seg| -> bool {
                let c0 = line_seg[0];
                let c1 = line_seg[1];
                if c0.y < c1.y {
                    p.y >= c0.y && p.y < c1.y
                } else if c0.y > c1.y {
                    p.y >= c1.y && p.y < c0.y
                } else {
                    false
                }
            }) {
                let c0 = line_seg[0];
                let c1 = line_seg[1];
                let a = c0.x - p.x + (p.y - c0.y) * (c1.x - c0.x) / (c1.y - c0.y);
                if (a > 0.0) {
                    println!("p is left of the line segment");
                } else {
                    println!("p is to the right of the line segment");
                }

                println!("{:?} {:?} {:?}", p, c0, c1);
                // now that the candidates are filtered, calculate intersections
        }
        false
    }

    fn gen_area(border: &Vec<Vert>) -> Vec<Tri> {
        let mut area = Vec::<Tri>::new();

        let triangles = triangulate(&Vert::verts_2_points(border.clone())).triangles;
        area.reserve(triangles.len() / 3);
        for (c0, c1, c2) in triangles.into_iter().tuples::<(_, _, _)>() {
            let mid = Vert {
                x: (border[c0].x + border[c1].x + border[c2].x) / 3.0,
                y: (border[c0].y + border[c1].y + border[c2].y) / 3.0
            };
            Domain::is_inside(&border, &mid);

            area.push(Tri::new(border, [c0, c1, c2]));
        }

        area
    }

    #[allow(dead_code, unused_variables)]
    fn new_from_rect(lower: Vert, upper: Vert) -> Domain {
        let mut border: Vec<Vert> = Vec::new();
        border.push(lower.clone());
        border.push(Vert {
            x: lower.x,
            y: upper.y,
        });
        border.push(upper.clone());
        border.push(Vert {
            x: upper.x,
            y: lower.y,
        });

        let mut domain = Domain {
            border,
            verts: Vec::<Vert>::new(),
            area: Vec::<Tri>::new(),
            cells: Vec::<SimCell>::new(),
        };

        domain.area.clone_from(&Domain::gen_area(&domain.border));

        domain
    }

    fn new_from_border(border: Vec<Vert>) -> Domain {
        let mut domain = Domain {
            border,
            verts: Vec::<Vert>::new(),
            area: Vec::<Tri>::new(),
            cells: Vec::<SimCell>::new(),
        };

        domain.area.clone_from(&Domain::gen_area(&domain.border));

        domain
    }

    #[allow(dead_code, unused_variables)]
    fn to_vtk(&self) -> Result<(), std::io::Error> {
        let mut file: File = File::create("domain.vtu")?;
        // Header
        file.write(
            format!(
                "<?xml version=\"1.0\"?>\n\
                   <VTKFile type=\"UnstructuredGrid\" version=\"0.1\" \
                   byte_order=\"LittleEndian\">\n  \
                     <UnstructuredGrid>\n    \
                         <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n    \
                         <Points>\n      \
                         <DataArray type=\"Float64\" Name=\"Points\" \
                         NumberOfComponents=\"3\" format=\"ascii\">\n        ",
                self.border.len(),
                self.area.len()
            )
            .as_bytes(),
        )?;
        // Point array
        for point in self.border.iter() {
            file.write(format!("{} {} {} ", point.x, point.y, 0.0f64).as_bytes())?;
        }

        // Point connectivity for cells
        file.write(
            b"\n      \
            </DataArray>\n    \
            </Points>\n    \
            <Cells>\n      \
            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n        ",
        )?;
        for tri in self.area.iter() {
            file.write(format!("{} {} {} ", tri.0, tri.1, tri.2).as_bytes())?;
        }

        // Offsets to ends of Cells
        file.write(
            b"\n      \
            </DataArray>\n      \
            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n        ",
        )?;
        for i in 0..self.area.len() {
            file.write(format!("{} ", (i + 1) * 3).as_bytes())?;
        }

        // Types
        file.write(
            b"\n      \
            </DataArray>\n      \
            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n        ",
        )?;
        // Triangles (Type 5)
        for _ in 0..self.area.len() {
            file.write(b"5 ")?;
        }
        file.write(
            b"\n      \
                </DataArray>\n    \
            </Cells>\n",
        )?;

        file.write(
            b"    \
            <CellData>\n      \
            <DataArray type=\"UInt32\" format=\"ascii\" \
                Name=\"AreaIndex\">\n        ",
        )?;
        for i in 0..self.area.len() {
            file.write(format!("{} ", i).as_bytes())?;
        }
        file.write(
            b"\n      \
            </DataArray>\n    \
            </CellData>\n",
        )?;

        file.write(
            b"\
            </Piece>\n\
            </UnstructuredGrid>\n\
            </VTKFile>\n",
        )?;

        Ok(())
    }
}

fn main() {
    let border = vec![
        Vert {x:  -9.30, y:  7.58 },
        Vert {x:   6.72, y:  7.62 },
        Vert {x:   7.00, y: -1.00 },
        Vert {x:  -1.58, y: -7.18 },
        Vert {x: -10.50, y: -4.18 },
        Vert {x:  -5.14, y: -1.24 },
        Vert {x:  -3.90, y:  2.82 },
    ];
    let domain = Domain::new_from_border(border);

    // Right now, this triangulates the points.
    // I want to remove the tris outside of the border
    domain.to_vtk().unwrap();
}
