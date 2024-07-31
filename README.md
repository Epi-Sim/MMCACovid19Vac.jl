# EpiSim.jl
A julia-based simulator for epidemic processes

## Using EpiSim.py Interface

The `EpiSim.py` file provides a Python interface for the EpiSim.jl epidemic simulator. Here's how to use it:

1. Import the `MMCACovid19` class:
   ```python
   from py_interface.EpiSim import MMCACovid19
   ```

2. Initialize the model:
   ```python
   model = MMCACovid19(
       executable_path="path/to/episim",
       config="path/to/config.json",
       data_folder="path/to/data",
       instance_folder="path/to/runs",
       initial_conditions="path/to/initial_conditions.nc"
   )
   ```

3. Run a single step:
   ```python
   new_state, next_date = model.step(start_date="2023-01-01", length_days=5)
   ```

4. Update configuration:
   ```python
   new_config = {...}  # Dict or path to JSON file
   model.update_config(new_config)
   ```

5. Run multiple steps:
   ```python
   current_date = "2023-01-01"
   for _ in range(10):
       new_state, next_date = model.step(start_date=current_date, length_days=5)
       current_date = next_date
   ```

Each model instance has a unique UUID and folder for storing state, allowing multiple simultaneous runs.

For more detailed examples, refer to the `run_model_example()` and `agent_flow_example()` functions in the `EpiSim.py` file.