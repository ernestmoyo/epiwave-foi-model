
#!/bin/bash
# ============================================
# Script to pull Nick Golding's ODE solver branch
# ============================================

cd "C:/Users/ernes/Documents/R/epiwave-foi-model"

echo "=== Adding Nick's fork as remote ==="
git remote add goldingn https://github.com/goldingn/epiwave-foi-model.git 2>/dev/null || echo "Remote already exists"

echo ""
echo "=== Fetching Nick's branch ==="
git fetch goldingn test-ode-solvers

echo ""
echo "=== Creating local branch to review ==="
git checkout -b review-nicks-ode-solvers goldingn/test-ode-solvers

echo ""
echo "=== Done! You are now on Nick's branch ==="
git branch
git log --oneline -10

echo ""
echo "=== Files in the repo ==="
ls -la

