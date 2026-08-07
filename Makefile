.PHONY: verify verify-python verify-integrity verify-lean paper main-ready clean

verify: verify-python verify-integrity

verify-python:
	python3 verify_exact.py
	python3 rank_packet_obstruction_tests.py
	python3 research/riemann_roch_visibility_certificate.py
	python3 research/torsion_packet_capacity_certificate.py
	python3 research/cyclic_trisection_trace_code_certificate.py
	python3 research/j0_six_channel_capacity_certificate.py
	python3 research/j0_cubic_family_rank_ceiling_certificate.py
	python3 research/j0_minimal_scale_obstruction_certificate.py
	python3 -m unittest discover -s tests -v
	git diff --exit-code

verify-integrity:
	test -s MANIFEST.sha256
	python3 scripts/build_integrity_manifest.py --check MANIFEST.sha256

verify-lean:
	test -s lake-manifest.json
	lake build

paper:
	cp paper/source/main.tex paper/paper.tex
	rm -rf paper/build
	mkdir -p paper/build
	latexmk -pdf -interaction=nonstopmode -halt-on-error \
		-output-directory=paper/build paper/paper.tex
	! grep -E 'Overfull \\hbox|Overfull \\vbox' paper/build/paper.log
	! grep -E 'LaTeX Warning:.*undefined' paper/build/paper.log
	cp paper/build/paper.pdf paper/paper.pdf
	rm -rf paper/build

main-ready: verify verify-lean paper

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf paper/build
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls \
		paper/*.fdb_latexmk paper/*.synctex.gz
