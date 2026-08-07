.PHONY: verify paper clean

verify:
	python3 verify_exact.py
	python3 rank_packet_obstruction_tests.py
	python3 research/riemann_roch_visibility_certificate.py
	python3 research/torsion_packet_capacity_certificate.py
	python3 -m unittest discover -s tests -v
	sha256sum -c MANIFEST.sha256

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error paper.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error paper.tex

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk paper/*.synctex.gz
