.PHONY: help install install-dev test lint format clean build
.DEFAULT_GOAL := help

help: ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install package in current environment
	pip install .

install-dev: ## Install package with development dependencies
	pip install -e ".[dev]"

test: ## Run tests
	pytest tests/ -v --cov=three_s_two_b --cov-report=term-missing

test-fast: ## Run tests without coverage
	pytest tests/ -v

lint: ## Run linting checks
	ruff check src/
	black --check src/
	mypy src/

format: ## Format code
	black src/ tests/
	ruff check src/ --fix

clean: ## Clean build artifacts
	rm -rf build/ dist/ *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean ## Build package
	python -m build

demo: ## Run a demo workflow
	@echo "Running 3S2B demo..."
	3s2b fragment "CCO" --max-mult 2
	@echo "Demo completed!"