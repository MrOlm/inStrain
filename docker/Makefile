NAME=instrain
# TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
# TAG="MO_dev"
VERSION:=$(shell cat ../inStrain/_version.py | awk -F'=' '{print $$2}')
VERSION:=$(shell echo $(VERSION)) # Gets rid of quotes

REGISTRY=mattolm/$(NAME)

all: build push run

build:
	docker image build --no-cache -t $(NAME):$(VERSION) -t $(NAME):latest -f Dockerfile .
	# yes | cp -a ../common_utils/ ./

	docker image build -t $(NAME):$(VERSION) -t $(NAME):latest -f Dockerfile .
	docker image tag $(NAME):latest $(REGISTRY):latest
	docker image tag $(NAME):$(VERSION) $(REGISTRY):$(VERSION)

push:
	docker image push $(REGISTRY):$(VERSION)
	docker image push $(REGISTRY):latest
	echo "mattolm/$(NAME):$(VERSION)" > LATEST

run:
	docker container run --rm $(NAME):latest
