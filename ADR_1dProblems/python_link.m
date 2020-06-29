% Useful links describing how to correctly load Python environments and
% Python modules in Matlab
% https://github.com/altermarkive/Calling-Python-from-Matlab
% https://www.mathworks.com/help/matlab/ref/pyenv.html
% https://www.mathworks.com/matlabcentral/answers/466974-how-to-change-python-path
% https://www.mathworks.com/help/matlab/matlab_external/undefined-variable-py-or-function-py-command.html

pyExec	= '/opt/anaconda3/bin/python3';
[ver, exec, loaded]	= pyversion(pyExec);

python_path = py.sys.path;

% Add folders to python system path to loead version 3.7 of Python.
if count(python_path, pyExec) == 0
    insert(py.sys.path, int64(0), pyExec);
end

% Verify that the Python 3.7 version from Anaconda is successfully loaded
[ver2, exec2, loaded2]	= pyversion;
assert(loaded2==1);

python_module = py.importlib.import_module('model_load_script');
python_module2 = py.importlib.import_module('hello');

inserted = zeros(2,2);
inserted(1,1) = 1;
inserted(1,2) = 2;
inserted(2,1) = 3;
inserted(2,2) = 4;
released = python_module2.say_hi(inserted);
python_module.model_load()




