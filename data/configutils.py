def state_atom_constraints(pocket, state):
    key = "atomtoi" if state == "inactive" else "atomtoa"
    values = pocket.get(key, [])
    if isinstance(values, str):
        return [values]
    if not isinstance(values, (list, tuple)) or not all(isinstance(v, str) for v in values):
        raise ValueError(f"Invalid {key} constraints")
    return list(values)
