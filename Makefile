.PHONY: verify-fast verify-all audit-deps audit-n64-kkt clean

verify-fast:
	bash ./scripts/run_fast.sh

verify-all:
	bash ./scripts/run_all.sh

audit-deps:
	python3 -m pip install -r requirements-audit.txt

audit-n64-kkt:
	python3 audits/n64/kkt_high_precision.py

clean:
	rm -f cases/n64/n64_code_fixed
