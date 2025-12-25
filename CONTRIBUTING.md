# Contributing to Charge Optimizer

Thank you for your interest in contributing! This project welcomes contributions from everyone.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your system info (OS, compiler version)

### Suggesting Features

Feature suggestions are welcome! Please open an issue describing:
- The feature you'd like to see
- Why it would be useful
- How it might work

### Code Contributions

1. **Fork the repository**

2. **Create a feature branch**
   ```bash
   git checkout -b feature/my-new-feature
   ```

3. **Make your changes**
   - Follow the existing code style
   - Add comments for complex logic
   - Update documentation if needed

4. **Test your changes**
   ```bash
   cd build
   make
   ctest
   ```

5. **Commit with clear messages**
   ```bash
   git commit -m "Add feature: brief description"
   ```

6. **Push and create a pull request**

## Code Style

- C++17 standard
- Use meaningful variable names
- Comment complex algorithms
- Keep functions focused and small
- Use Eigen for linear algebra

## Areas for Contribution

### High Priority
- [ ] Interior-point QP solver
- [ ] ADMM solver for large problems
- [ ] Python bindings (pybind11)
- [ ] Multi-conformer fitting
- [ ] GPU acceleration (CUDA)

### Medium Priority
- [ ] Additional force field support
- [ ] Web interface
- [ ] More validation metrics
- [ ] Performance benchmarks
- [ ] Integration with molecular visualization tools

### Good First Issues
- [ ] Add more example molecules
- [ ] Improve error messages
- [ ] Add unit tests
- [ ] Documentation improvements
- [ ] Code cleanup and refactoring

## Testing

All contributions should include appropriate tests. Run the test suite:

```bash
cd build
ctest --output-on-failure
```

## Documentation

Update the README.md if you:
- Add new features
- Change the API
- Add new dependencies

## Questions?

Feel free to open an issue for any questions about contributing!
