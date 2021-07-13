use std;
extern crate gurobi;
use gurobi::*;


fn main() {  

    let env = Env::new("logfile.log").unwrap();
    let mut model = env.new_model("model1").unwrap();
    
    let mut x = model.add_var("x", Continuous, 2.0, 0.0, INFINITY, &[], &[]).unwrap();
    let mut y = model.add_var("y", Continuous, 3.0, 0.0, INFINITY, &[], &[]).unwrap();

    let mut obj_expression: LinExpr = LinExpr::new();
    obj_expression = obj_expression.add_term( 2.0, x.clone() );
    obj_expression = obj_expression.add_term( 3.0, y.clone() );
    let mut constraint = LinExpr::new();
    constraint = constraint.add_term(1.0, x.clone());
    constraint = constraint.add_term(1.0, y.clone());

    model.add_constr("constraint", constraint, Equal, 1.0);
    model.add_constr("constraint", x.clone() + 0.5, Greater, 1.0);
    model.update().unwrap();
    model.set_objective(obj_expression,Maximize).unwrap();
    model.write("logfile_new.lp").unwrap();
    model.optimize().unwrap(); // Solve result
    let res = model.get_values(attr::X, &[x,y]).unwrap();
    

    println!("x :{} y:{}",res[0],res[1]);

}