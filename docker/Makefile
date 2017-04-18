IMAGE=ont-assembly-polish

build: Dockerfile
	docker build --no-cache --pull --build-arg CACHEBUST=$(date +%s) -t $(IMAGE) .
