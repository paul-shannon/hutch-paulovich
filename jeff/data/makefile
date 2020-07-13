ID = 13jul2020
APP_IMAGE_NAME = jeff_app_$(ID)

build:
	docker build -t $(APP_IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -p 3838:3838 \
           --entrypoint /bin/bash \
           $(APP_IMAGE_NAME)

