use std;
extern crate gurobi;
use gurobi::*;

// This file illustrate how to set up a model with the style of the example on the gurobi crate website

fn main() {  
    let env = Env::new("logfile.log").unwrap();
    let mut model = env.new_model("model1").unwrap();

    // initialize the vector: x+
    let mut x_pos = Vec::new();
        
    for i in 0..7{

        let  name = format!("{}{}", "x_pos", i);

        let  str_name = &name[..];
        // Should put weight in objective function because we don't want weight in constraint
        x_pos.push(model.add_var(str_name, Continuous,1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: x-
    let mut x_neg = Vec::new();
    
    for i in 0..7{

        let  name = format!("{}{}", "x_neg", i);

        let  str_name = &name[..];
        x_neg.push(model.add_var(str_name, Continuous, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }
    let q0 = model.add_var("q0", Continuous , 0.0, 0.0, INFINITY,  &[],  &[]).unwrap();
    let q1 = model.add_var("q1", Continuous , 0.0, 0.0, INFINITY,  &[],  &[]).unwrap();

    model.set_objective(2.0 * &x_pos[0]+ 2.0 * &x_pos[1]+ 2.0 * &x_pos[2]
        + &x_pos[3]+ 2.0 * &x_pos[4]+ &x_pos[5]+ &x_pos[6] + 2.0 * &x_neg[0]+ 2.0 * &x_neg[1]+ 2.0 * &x_neg[2]+ &x_neg[3] 
        + 2.0 * &x_neg[4]+ &x_neg[5]+ &x_neg[6], Minimize);
    
    model.add_constr("constraint 1",    &x_pos[0] - &x_neg[0] , Equal, -1.0);
    model.add_constr("constraint 2",    &x_pos[1] - &x_neg[1] , Equal, 1.0);
    model.add_constr("constraint 3",    &x_pos[2] - &x_neg[2] - 1.0 * &q1 , Equal, -1.0);
    model.add_constr("constraint 4",    &x_pos[3] - &x_neg[3] + 1.0 * &q1, Equal, 0.0);
    model.add_constr("constraint 5",    &x_pos[4] - &x_neg[4] - 1.0 * &q0, Equal, 1.0);
    model.add_constr("constraint 6",    &x_pos[5] - &x_neg[5] + 1.0 * &q0, Equal, 0.0);
    model.add_constr("constraint 7",    &x_pos[6] - &x_neg[6] - 1.0 * &q0 - 1.0 * &q1, Equal, 0.0);
    model.update().unwrap();
    model.write("logfile_2.lp").unwrap();
    model.optimize().unwrap(); // Solve result

    let x_neg_val = model.get_values(attr::X, &x_neg).unwrap(); // Get the result for x+ and x-
    let x_pos_val = model.get_values(attr::X, &x_pos).unwrap();

    println!("x_neg_val :{:?}",x_neg_val);
    println!("x_pos_val :{:?}",x_pos_val);



}








