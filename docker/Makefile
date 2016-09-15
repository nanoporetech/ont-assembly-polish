IMAGE=ont-assembly-polish

build: Dockerfile
	docker build --no-cache --pull --build-arg CACHEBUST=$(date +%s) -t $(IMAGE) .

.PHONY: push
push:
	docker tag $(IMAGE):latest docker-registry.oxfordnanolabs.local:5000/$(IMAGE):latest
	docker push docker-registry.oxfordnanolabs.local:5000/$(IMAGE)
