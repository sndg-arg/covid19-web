SNDG-COVID19
============

SNDG web to show argentinian COVID19 isolates

# Build from the scratch

## Preconditions for the tutorial
    * Python 3.7 + Pip , VirtualEnv recommended
    * Docker
    * Environment variables
        * export DJANGO_READ_DOT_ENV_FILE=True
        * $DB_PASS
        * .env
        
## Web
### Install JBrowse

```console
    #https://jbrowse.org/docs/embedding.html
    wget "https://github.com/GMOD/jbrowse/releases/download/1.16.6-release/JBrowse-1.16.6.zip"
    unzip JBrowse-1.16.6.zip
    rm JBrowse-1.16.6.zip
    mv JBrowse-1.16.6/ data/jbrowse
    sed -i 's|c.src=function(t){return i.p+""+|c.src=function(t){return "/static/jbrowse/"+i.p+""+|' data/jbrowse/dist/browser.bundle.js
```
    
### Build Front End libraries
```console
    # Install Javascript dependencies
    docker run -u $(id -u):$(id -g) --rm -v $PWD:/out -w /out node:lts npm install
    # Fix old libraries
    sed -i 's|require("bootstrap/js/tooltip.js")|require("bootstrap/js/dist/tooltip.js")|' ./node_modules/feature-viewer/lib/index.js
    sed -i 's|require("bootstrap/js/popover.js")|require("bootstrap/js/dist/popover.js")|' ./node_modules/feature-viewer/lib/index.js
    # Build javascript dependency
    docker run -u $(id -u):$(id -g) --rm -v $PWD:/out -w /out node:lts npm run-script build
```

* DB
    * Postgres Engine
        <code>
        docker run -v $PWD/data/sndgr:/var/lib/postgresql/data --name sndgr \
        -e POSTGRES_PASSWORD=$DB_PASS -e POSTGRES_DB=sndgr -p 5432:5432 -d \
        --shm-size 512m postgres:11
        </code>
    * Base data
        * ./manage.py migrate
        * ./manage.py init_db
        * ./manage.py import_genome -i data/curated/covid19.gb -a "COVID19" -n "COVID19" -t 2697049
        * ./manage.py update_covid_data

    * Annotation
        * ./manage.py pdbfiles2pdbdb # change PDBfiles path with "pdbs_dir" parameters
        * ./manage.py fpocket2pdbdb
    * Load genome in JBrowse
        * ./manage.py db2jbrowse --dbname COVID19
    * Create admin user
        * ./manage.py createsuperuser

* Start server
     ./manage.py runserver

* Development
    * ./manage.py makemessages  -l es -l en  -i -v 3 "env/*" # env in case you have an virtualenv there
    * manage.py shell_plus --ipython --print-sql # interactive ipython shell
