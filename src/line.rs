/*a Copyright

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@file    bezier_line.rs
@brief   Part of geometry library
 */

//a Imports
use geo_nd::{Float, Vector};
use crate::Bezier;

//a BezierLineIter
//tp BezierLineIter
/// An iterator of straight lines that form a single Bezier curve
///
/// An iteration will provide (Pa, Pb) pairs of points, with
/// the next iteration providing (Pb, Pc), then (Pc, Pd), etc;
/// sharing the end/start points.
pub struct BezierLineIter<F:Float, V:Vector<F,D>, const D:usize> {
    /// Maximum curviness of the line segments returned
    straightness: F,
    /// A stack of future beziers to examine
    /// The top of the stack is p0->p1; below that is p1->p2, etc
    /// These beziers will need to be split to achieve the maximum
    /// curviness
    stack : Vec<Bezier<F,V,D>>
}

//pi BezierLineIter
impl <F, V, const D:usize> BezierLineIter<F, V, D>
where F:Float, V:Vector<F,D> {
    //fp new
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn new(bezier:&Bezier<F, V, D>, straightness:F) -> Self {
        let mut stack = Vec::new();
        stack.push(bezier.clone());
        Self { straightness, stack }
    }

    //zz All done
}

//ip Iterator for BezierLineIter
impl <F, V, const D:usize> Iterator for BezierLineIter<F, V, D>
where F:Float, V:Vector<F,D> {
    /// Item is a pair of points that make a straight line
    type Item = (V, V);
    /// next - return None or Some(pa,pb)
    ///
    /// It pops the first Bezier from the stack: this is (pa,px); if
    /// this is straight enough then return it, else split it in two
    /// (pa,pm), (pm,px) and push them in reverse order, then recurse.
    ///
    /// This forces the segment returned (eventually!) to be (pa,pb)
    /// and to leave the top of the stack starting with pb.
    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop() {
            None => None,
            Some(b) => {
                if b.is_straight(self.straightness) {
                    Some(b.endpoints())
                } else {
                    let (b0, b1) = b.bisect();
                    self.stack.push(b1);
                    self.stack.push(b0);
                    self.next()
                }
            }
        }
    }

    //zz All done
}

