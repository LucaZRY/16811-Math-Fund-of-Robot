# HW1 LDU Package

This folder contains a clean, reusable setup for your HW1 LDU decomposition.

## Files
- `LDU_Decomposition.m` — Standalone MATLAB function implementing PA = L*D*U with partial pivoting.
- `HW1.m` — Script that demonstrates calling the function (convert to Live Script if desired).
- `test_hw1.m` — Quick unit tests across multiple matrices.

## How to Use
1. Open MATLAB and set your Current Folder to this directory.
2. Run `HW1.m` to see a demo.
3. (Optional) Convert `HW1.m` to a Live Script: **Home → Save As → Live Script (.mlx)**.
4. Run `test_hw1.m` for a quick verification across several matrices.

## Notes
- Keep `LDU_Decomposition.m` in the same folder as your script/live script, or add the folder to your MATLAB path.
- The function returns P, L, D, U such that **P*A = L*D*U**, with L unit-lower, U unit-upper (diagonal factored into D).
