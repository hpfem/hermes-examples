Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
    MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(prev_slns, KAPPA));
    MeshFunctionSharedPtr<double>  pressure(new PressureFilter(prev_slns, KAPPA));

    ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
    ScalarView Mach_number_view("Mach number", new WinGeom(650, 0, 600, 300));
    ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
    ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
    OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
  #pragma endregion

EulerEquationsWeakFormStabilization wf_stabilization(prev_rho);
DiscreteProblem<double> dp_stabilization(&wf_stabilization, space_stabilization);

if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
  wf.set_stabilization(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, NU_1, NU_2);

LinearSolver<double> solver(&wf, spaces);

#pragma region 5. Time stepping loop.
int iteration = 0;
for(double t = 0.0; t < TIME_INTERVAL_LENGTH; t += time_step_n)
{
  // Info.
  Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

  if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
  {
    int mesh_size = space_stabilization->get_num_dofs();
    assert(mesh_size == space_stabilization->get_mesh()->get_num_active_elements());
    dp_stabilization.assemble(rhs_stabilization);
    if(!wf.discreteIndicator)
    {
      wf.set_discreteIndicator(new bool[mesh_size], mesh_size);
      memset(wf.discreteIndicator, 0, mesh_size * sizeof(bool));
    }
    Element* e;
    for_all_active_elements(e, space_stabilization->get_mesh())
    {
      AsmList<double> al;
      space_stabilization->get_element_assembly_list(e, &al);
      if(rhs_stabilization->get(al.get_dof()[0]) >= 1)
        wf.discreteIndicator[e->id] = true;
    }
  }

  // Set the current time step.
  wf.set_current_time_step(time_step_n);

  try
  {
    // Solve.
    solver.solve();

#pragma region *. Get the solution with optional shock capturing.
    if(!SHOCK_CAPTURING)
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, prev_slns);
    else
    {
      FluxLimiter* flux_limiter;
      if(SHOCK_CAPTURING_TYPE == KUZMIN)
        flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), spaces, true);
      else
        flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), spaces);

      if(SHOCK_CAPTURING_TYPE == KUZMIN)
        flux_limiter->limit_second_orders_according_to_detector();

      flux_limiter->limit_according_to_detector();

      flux_limiter->get_limited_solutions(prev_slns);
    }
#pragma endregion
  }
  catch(std::exception& e) { std::cout << e.what(); }

  // Calculate time step according to CFL condition.
  CFL.calculate(prev_slns, mesh, time_step_n);

#pragma region 5.1. Visualization
  if((iteration - 1) % EVERY_NTH_STEP == 0) 
  {
    // Hermes visualization.
    if(HERMES_VISUALIZATION)
    {        
      Mach_number->reinit();
      pressure->reinit();
      pressure_view.show(pressure, 1);
      Mach_number_view.show(Mach_number, 1);
      order_view.show(space_rho);
    }
    // Output solution in VTK format.
    if(VTK_VISUALIZATION)
    {
      pressure->reinit();
      Linearizer lin;
      char filename[40];
      sprintf(filename, "Pressure-%i.vtk", iteration - 1);
      lin.save_solution_vtk(pressure, filename, "Pressure", false);
      sprintf(filename, "VelocityX-%i.vtk", iteration - 1);
      lin.save_solution_vtk(prev_rho_v_x, filename, "VelocityX", false);
      sprintf(filename, "VelocityY-%i.vtk", iteration - 1);
      lin.save_solution_vtk(prev_rho_v_y, filename, "VelocityY", false);
      sprintf(filename, "Rho-%i.vtk", iteration - 1);
      lin.save_solution_vtk(prev_rho, filename, "Rho", false);
    }
  }
#pragma endregion
}
#pragma endregion

// Done.
return 0;