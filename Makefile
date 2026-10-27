.PHONY: check test figures clean

PROJECT := projects/high-energy-water

check:
	$(MAKE) -C $(PROJECT) check

test:
	$(MAKE) -C $(PROJECT) test

figures:
	$(MAKE) -C $(PROJECT) figures

clean:
	$(MAKE) -C $(PROJECT) clean

