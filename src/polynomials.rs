
use crate::{DPLUS2_CHOOSE_2, algebraic_types::Matrix, DEGREE, FIELD_ORDER, COEFF_BIT_SIZE};

pub fn internal_add(a: u64,b: u64) -> u64 {
  match FIELD_ORDER {
    2 => internal_add_f2(a,b),
    3 => internal_add_f3_fast(a,b),
    _ => panic!("Field size not supported"),
  }
}

pub fn internal_add_f2(a: u64,b: u64) -> u64 {
  a ^ b
}

#[allow(dead_code)]
pub fn internal_add_f3(a: u64, b: u64) -> u64 {
  const M1: u64 = 0x5555555555555555;
  const M2: u64 = 0xAAAAAAAAAAAAAAAA;

  let xor = a^b;
  let and = a&b;

  let one = (and & M1) << 1;
  let two = (and & M2) >> 1;

  let ab = ((a&M2) >> 1) & b;
  let ba = ((b&M2) >> 1) & a;

  let mul = (ab | ba) * 0b11;

  (mul ^ xor) | one | two
}

pub fn internal_add_f3_fast(a: u64,b: u64) -> u64 {
  const M2: u64 = 0xAAAAAAAAAAAAAAAA; 
  let na=!a;
  let nb=!b;
  let a4= ((M2 & na) >> 1) & na;
  let b4= ((M2 & nb) >> 1) & nb;
  !((  (a4 << 1 | a4) | (b4 << 1 | b4))^(a|b))
} 


#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Polynomial {
  pub bits: u64
}

impl Polynomial {
  pub fn new(bits: u64) -> Polynomial {
    Polynomial { bits: bits }
  }

  pub fn str(&self, lut: &Vec<Term>) -> String {
    let mut poly_str = String::new();
    let mut empty = true;
    for i in 0..DPLUS2_CHOOSE_2 {
      const MASK: u64 = !(!0 << COEFF_BIT_SIZE); // only the last COEFF_BIT_SIZE bits are 1.
      let constant = ((self.bits >> (i*COEFF_BIT_SIZE)) & MASK) as u8;
      if constant > 0 {
        if empty {
          poly_str = format!("{}", lut[i].create_similar(constant).str());
          empty = false;
        } else {
          poly_str = format!("{} {}", poly_str, lut[i].create_similar(constant).str());
        }
      }
    }
    poly_str
  }

  #[allow(dead_code)]
  pub fn print(&self, lut: &Vec<Term>) {
    println!("{}", self.str(lut));
  }

  pub fn generate_default_lut() -> Vec<Term> {
    (0..=DEGREE as u8)
        .flat_map(move |a| {
            (0..=(DEGREE as u8) - a).map(move |b| {
                let c = DEGREE as u8 - b - a;
                Term {
                    x_deg: a,
                    y_deg: b,
                    z_deg: c,
                    constant: 1,
                }
            })
        })
        .collect()
  }

  pub fn transform_by_matrix(self, transform_lut: &Vec<u64>) -> Polynomial {
    let mut bits = 0;
    for i in 0..DPLUS2_CHOOSE_2 {
      let coeff = (self.bits >> (i*COEFF_BIT_SIZE)) & (!(!0 << COEFF_BIT_SIZE));
      if coeff > 0 {
        let new_bits = Polynomial::multiply_bits_by_constant(transform_lut[i], coeff);
        bits = internal_add(bits, new_bits);
      }
    }
    Polynomial { bits: bits }
  }

  pub fn multiply_bits_by_constant(bits: u64, constant: u64) -> u64 {
    match constant % FIELD_ORDER as u64 {
      0 => 0,
      1 => bits,
      _ => {
        match FIELD_ORDER {
          2 => bits,
          3 => Polynomial::multiply_bits_by_2_f3(bits),
          _ => panic!("Field size not supported"),
        }
      }
    }
  }

  fn is_leading_coeff_unit(poly: Polynomial) -> bool {
    match FIELD_ORDER {
      2 => true,
      3 => poly.bits.leading_zeros() % 2 == 1,
      _ => panic!("Field size not supported"),
    }
  }
  
  pub fn multiply_leading_coeff_to_unit(self) -> Polynomial {
    match FIELD_ORDER {
      2 => self,
      3 => {
        if !Polynomial::is_leading_coeff_unit(self) {
          Polynomial::new(Polynomial::multiply_bits_by_2_f3(self.bits))
        } else {
          self
        }
      },
      _ => panic!("Field size not supported"),
    }
  }

  pub fn mul_constant(self, constant: u64) -> Polynomial {
    Polynomial { bits: Polynomial::multiply_bits_by_constant(self.bits, constant) }
  }
    

  pub fn multiply_bits_by_2_f3(bits: u64) -> u64 {
    const M1: u64 = 0x5555555555555555; 
    const M2: u64 = 0xAAAAAAAAAAAAAAAA;
    let a1 = (bits & M2) >> 1;
    let a2 = (bits & M1) << 1;
    a1 | a2
  }
}



#[derive(Debug, Copy, Clone, PartialEq, Hash, Eq)]
pub struct Term { 
  pub x_deg: u8,
  pub y_deg: u8,
  pub z_deg: u8,
  pub constant: u8,
}


impl Term {
  pub fn str(self) -> String {
    format!("{}_{}{}{}", self.constant, self.x_deg, self.y_deg, self.z_deg)
  }

  pub fn transform_by_matrix(self, matrix: &Matrix, lut: &Vec<Term>) -> u64 {
    if self.constant == 0 {
      return 0;
    }
    let p1 = exponentiate_linear_polynomial(matrix.data[0][0], matrix.data[0][1], matrix.data[0][2], self.x_deg);
    // print_vec_term(&p1);
    let p2 = exponentiate_linear_polynomial(matrix.data[1][0], matrix.data[1][1], matrix.data[1][2], self.y_deg);
    // print_vec_term(&p2);
    let p3 = exponentiate_linear_polynomial(matrix.data[2][0], matrix.data[2][1], matrix.data[2][2], self.z_deg);
    // print_vec_term(&p3);
    let terms = polynomial_product(p1, p2, p3);
    // print_vec_term(&terms);
    // TESTING CONLUSION: term degrees are good, unsure of constant.

    let mut result = 0;
    for t in terms {
      for i in 0..DPLUS2_CHOOSE_2 {
        if lut[i].is_similar(t) { // TODO: Not this idiotic way of finding the right term index.
          let bits = (t.constant as u64) << (i*COEFF_BIT_SIZE);
          result = internal_add(result, bits);
        }
      }
    }
    result
  }

  pub fn is_similar(&self, t: Term) -> bool {
    self.x_deg == t.x_deg && self.y_deg == t.y_deg && self.z_deg == t.z_deg
  }

  pub fn create_similar(&self, constant: u8) -> Term {
    Term { x_deg: self.x_deg, y_deg: self.y_deg, z_deg: self.z_deg, constant: constant }
  }
}

#[allow(dead_code)]
pub fn print_vec_term(v: &Vec<Term>) {
  for t in v {
    print!("{} +", t.str());
  }
  println!();
}

// (ax+by+cz)^m = sum_{k1+k2+k3=m} (m choose k1,k2,k3) a^k1 b^k2 c^k3 x^k1 y^k2 z^k3
// SLOW!
pub fn exponentiate_linear_polynomial(a: u8, b: u8, c: u8, m: u8) -> Vec<Term> {
  let mut terms: Vec<Term> = Vec::new();
  for k1 in 0..=m {
    for k2 in 0..=(m-k1) {
      let k3 = m-k1-k2;
      let coeff = binomial_coefficient(m, k1, k2, k3) % FIELD_ORDER as u8;
      if (coeff == 0)|| (k1>0 && a==0) || (k2>0 && b==0) || (k3>0 && c==0) {
        continue;
      }
      let final_coeff = (coeff * a.pow(k1 as u32) * b.pow(k2 as u32) * c.pow(k3 as u32)) % FIELD_ORDER as u8;
      let term = Term { x_deg: k1, y_deg: k2, z_deg: k3, constant: final_coeff };
      terms.push(term);
    }
  }
  terms
}

pub fn polynomial_product(a: Vec<Term>, b: Vec<Term>, c: Vec<Term>) -> Vec<Term> {
  let mut result: Vec<Term> = Vec::new();
  for t1 in &a {
    for t2 in &b {
      for t3 in &c {
        let term = Term { x_deg: t1.x_deg + t2.x_deg + t3.x_deg, 
                          y_deg: t1.y_deg + t2.y_deg + t3.y_deg, 
                          z_deg: t1.z_deg + t2.z_deg + t3.z_deg, 
                          constant: (t1.constant * t2.constant * t3.constant) % FIELD_ORDER as u8}; // TODO: F4 wont work like this
        result.push(term);
      }
    }
  }
  result
}

pub fn binomial_coefficient(m: u8, k1: u8, k2: u8, k3:u8) -> u8 {
  (factorial(m) / (factorial(k1) * factorial(k2) * factorial(k3))) as u8
}

pub fn factorial(n: u8) -> u64 {
  let mut result: u64 = 1;
  for i in 1..=n as u64 {
    result *= i;
  }
  result
}

pub fn generate_transform_lut(pgl3: &Vec<Matrix>, lut: &Vec<Term>) -> Vec<Vec<u64>> {
  let mut result = vec![];
  for m in pgl3 {
    let mut result_for_m = vec![];
    for t in lut{
      let transformed = t.transform_by_matrix(&m, &lut);
      result_for_m.push(transformed);
    }
    result.push(result_for_m);
  }
  result
}