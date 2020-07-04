unit:
	pytest -v

coverage:
	pytest --cov=helperlibs --cov-report=html --cov-report=term-missing

lint:
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 . --count --exit-zero --max-complexity=20 --statistics
