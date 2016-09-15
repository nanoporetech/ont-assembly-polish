# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Create working directory if necessary:
WDT=$(WDIR)/PIPELINE_NAME # Depending on a file instead of the working directory to avoid spurious re-runs.
wdir: $(WDT)
$(WDT):
	@mkdir -p $(WDIR); echo $(PIPELINE_NAME) > $(WDT)

# Delete working directory:
.PHONY: clean_wdir
clean_wdir:
	@rm -fr $(WDIR)

# Delete results:
.PHONY: clean_res
clean_res:
	@rm -r $(RES)/

# Print pipeline info:
.PHONY: info
info:
	@echo Pipeline name: $(PIPELINE_NAME)
	@echo Pipeline working directory: $(WDIR)
	@echo Pipeline repository: $(REPO)

# Commit all changes:
.PHONY: com
com:
	@git commit -a

# Install python requirements:
.PHONY: dep
dep:
	pip install -r requirements.txt

# Print help:
.PHONY: help
help:
	@echo "Useful targets:"
	@echo "wdir			create working directory"
	@echo "clean_wdir		delete working directory. WARNING: all data will be lost!"
	@echo "clean_res		delete results directory. WARNING: all data will be lost!"
	@echo "info			print pipeline info"
	@echo "com 			commit all changes"
	@echo "dep 			install python dependencies using pip"
