# Thesis-Code

cd ~/Desktop/THESIS\ CODE

git status

git add name_of_code.py name_of_code.py

git add .

git commit -m "Explanation of changes"

git push

If too big then:

find . -type f -size +40M

git config --global http.postBuffer 524288000

echo "# test" > test_push.txt
git add test_push.txt
git commit -m "Test pushing a small file"
git push


# Quetzalcoatl

ssh giacomovargiu@quetzalcoatl.upc.edu

vPVrJ2DAvPFa

cd ~/thesis-code/Thesis-Code

git pull origin main


---

## Campaign: TOY EXAMPLE
Just a test to check if everything executes correctly.

### Launch
- Login into Quetz computer (only if trying to execute this on the HPC cluster).
- Move to repo directory: `cd ~/Thesis-Code/` 
- Activate virtual environment: `source venv/bin/activate`
- Launch main script: `python3 main_toy.py`

#### Advanced execution:
- If unattended execution is desired, use nohup command: `nohup python3 main_toy.py &`.
- If redirection to a log file is desired: `python3 main_toy.py > main_toy.log 2>&1`

Here is a one-liner that executes the campaign in headless mode redirecting the logs:
```shell 
nohup bash -lc 'cd ~/Thesis-Code && source venv/bin/activate && python3 main_toy.py' > ~/Thesis-Code/main_toy.log 2>&1 &
```

After running the command you’ll get a job number and PID, e.g.: `[1] 12345`.
You can verify it’s still running with: `ps -p 12345`.
To stop the job, use `kill 12345`, or if the PID is unknown use `pkill -f "python3 main_toy.py"`.

While the program is running, you can check what it is doing by following the log file:
```shell 
tail -f ~/Thesis-Code/main_toy.log
```

### Post-process
The campaign will generate several outputs inside the folder `output/toy_example/`.

---
