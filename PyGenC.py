
import sys
import subprocess
import os
import re
import shutil
import random
import threading
import tempfile
import regex

memprfx = ""
memsfx = "_"
maxexonlength = 100
sex = False
charlist = ' -+&\\/|!=()\"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
mutation_rate = 0.02
max_noncompile_gen = 3

class execution(object):
    def __init__(self, name):
        self.intron = None
        self.exon = None
        self.includes = []

# represents a gene function
class exon(object):

    def __init__(self, name):
        self.name = name
        self.introns = []
        self.expression = ""
        self.testdef = ""
        self.get = None # defines exon direct intron get
        self.defaultget = ""
        self.set = [] # defines exon direct intron sets
        self.public = False
        self.definition = ""
        self.executions = []
        self.gene = ""

    def resolve_executions(self, depth):

        self.definition = ""
        includes = []

        for e in self.executions:
            self.definition += depth+e.exon.name + " ("

            count = 0
            for i in e.exon.introns:
                self.definition += i.name
                if count != len(e.exon.introns)-1:
                    self.definition += i.seperator
                count+=1    
            self.definition += ")"

            # check for sub exe
            e.exon.resolve_executions(depth + "  ")
            

            if len(e.exon.executions) > 0:
                self.definition += depth + "{\n"
                self.definition += e.exon.definition
                self.definition += depth + "}\n"
            else:
                self.definition += ";\n"

            includes.extend(e.includes)
        return includes

class intron(object):

    def __init__(self, name, type, public):
        self.name = name
        self.type = type
        self.public = public
        self.defaultval = ""
        self.seperator = ", "

class gene(object):

    def __init__(self, name):
        self.name = name
        self.namespaces = []
        self.inheritance = []
        self.exons = []
        self.introns = []
        self.includes = []
        self.constdef = []

    def resolve_exons(self):
        # build execution string for each exon
        for e in self.exons:
             
             # reset execution
            self.includes.extend(e.resolve_executions("  "))

        self.includes = list(dict.fromkeys(self.includes))

        # check if self includes if so remove
        if self.name in self.includes:
            self.includes.remove(self.name)

        if "" in self.includes:
            self.includes.remove("")


def get_fitness(g, line):

    score = 0
    for e in g.exons:

        if e.name + "(" in line:
            score += 1
        if e.name in line:
            score += 1
        if e.expression in line:
            score += 1

    for i in g.introns:
        if i.name+" = " in line:
            score += 1
        if i.name in line:
            score += 1
        if i.type in line:
            score += 1
            
    return score

def addexon(g, line):

    result = re.search(',(.*)\)', line)
    cdef =  result.groups(1)
   
    if "Constructor" in cdef[0]:
        e = exon(g.name)
        if "CopyConstructor" in cdef[0]:
            e.introns.append(intron("rhs", g.name+"& ", False))
    else:
        e = exon(cdef[0].strip())
        e.expression = "void"

    g.exons.append(e)

def construct_repo_dir(name):

    #check to construct base_repo
    if not os.path.exists(os.path.dirname("./"+name+"/")):
        os.makedirs(os.path.dirname("./"+name+"/"))
    
    #check to construct source
    if not os.path.exists(os.path.dirname("./"+name+"/src/")):
        os.makedirs(os.path.dirname("./"+name+"/src/"))
    
    #check to construct include
    if not os.path.exists(os.path.dirname("./"+name+"/includes/")):
        os.makedirs(os.path.dirname("./"+name+"/includes/"))
    
    #check to construct tests
    if not os.path.exists(os.path.dirname("./"+name+"/tests/")):
        os.makedirs(os.path.dirname("./"+name+"/tests/"))

def create_testmain(genes, path):

    print("finalizing unittests... ")
    CodeString=""
    for g in genes:
        # include the unit test files
        CodeString += "#include \"../includes/" +g.name+".h\"\n"

    #include gtest header
    CodeString += "#include <gtest/gtest.h>\n\n"
    CodeString += "int main(int argc, char **argv)\n"
    CodeString += "{\ntesting::InitGoogleTest(&argc, argv);\n"
    CodeString += "return RUN_ALL_TESTS();\n}\n"

    cfile = open(path+"/tests/test_main.cpp", 'w')
    cfile.write(CodeString)
    cfile.close()


def create_cmakelists(genes, path):

    print("generating cmakelists... ")
    CodeString = "cmake_minimum_required(VERSION 2.6)\n\n"
    CodeString += "set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} -std=c++11\")\n"
    CodeString += "# Locate GTest\n"
    CodeString += "find_package(GTest REQUIRED)\n"
    CodeString += "include_directories(${GTEST_INCLUDE_DIRS})\n\n"
    CodeString += "# Link runTests with what we want to test and the GTest and pthread library\n"
    CodeString += "add_executable(executeTests tests/test_main.cpp"

    for g in genes:
        CodeString += " src/"+g.name+".cpp tests/"+g.name+"_test.cpp"

    CodeString += ")\ntarget_link_libraries(executeTests ${GTEST_LIBRARIES} pthread)"
    makefile = open(path+"/CMakeLists.txt", 'w')
    makefile.write(CodeString)
    makefile.close()

# rules
# 1) Each class with 1 header and cpp and 1 associated test
def create_test(g, path):
    print("generating test file... ", g.name+"_test" +(".cpp"))
    CodeString = "\n/**************** Auto Generated File **********************/\n\n"

    #include the header file
    CodeString += "#include \"../includes/" +g.name+".h\"\n"

    #include gtest header
    CodeString += "#include <gtest/gtest.h>\n\n"

    # create the test constdef
    for c in g.constdef:
        CodeString += c + "\n"

    # define tests for each exon
    for e in g.exons:
        if e.testdef != "":
            name = e.name

            # check default constructor and copy constructor
            if e.name == g.name and len(e.introns)==0:
                name = "DefaultConstructor"
            elif e.name == g.name and len(e.introns)==1 and g.name+"&" in e.introns[0].type:
                name = "CopyConstructor"
            elif e.name == g.name:
                name = "Constructor"
                for i in e.introns:
                    name += "_"+i.name

            CodeString += "TEST("+g.name+", "+name+") {\n"
            CodeString += e.testdef + "}\n\n"

    cfile = open(path+"/tests/"+ g.name + "_test.cpp", 'w')
    cfile.write(CodeString)
    cfile.close()

def gen_mutation(g):

    while True:
        line = ""
        
        charlength = random.randint(1, maxexonlength)
        index = 0

        while index < charlength:
            line += random.choice(charlist)
            index = index + 1
        line += ";"
        # try to filter invalid syntax here as much as possible
        if line.count("(") != line.count(")") or line[-1] != ";" or line.count('"') % 2 != 0 or line.count("{") != line.count("}") or ("=" in line and not " = " in line):
          #  print("Mutation Failed retrying ", line)
            continue
        else:
            print("Mutation Successful, inserting ", line)
            break

    return line

def getdefaultfromtype(type):
    
    if "*" in type:
        return "nullptr"
    if "int" in type:
        return "0"
    if "float" in type:
        return "0.0f"
    if "double" in type:
        return "0"

    return ""

def create_cpp(g, path):

    fitness = 0
    print("generating cpp... ", g.name +(".cpp"))
    CodeString = "\n/**************** Auto Generated File **********************/\n\n"

    #include header
    CodeString += "#include \"../includes/" +g.name+".h\"\n\n"
    
    # define constructors first
    for e in g.exons:
            
            if e.name != g.name:
              CodeString += e.expression + " "
            CodeString +=  g.name + "::" + e.name + "("
            for i in e.introns:
                CodeString += i.type + i.name
                if i != e.introns[-1]:
                    CodeString += i.seperator
            CodeString += ")"

            # check to initialize members with defaultget
            if e.name == g.name:

                if len(g.introns) > 0:
                    CodeString += ": \n"

                count = 0
                for i in g.introns:
                    type = getdefaultfromtype(i.type)

                    if type != "":
                        CodeString += i.name + "("+getdefaultfromtype(i.type)+")"
                        if count != len(g.introns)-1:
                            CodeString += "\n,"

                    count +=1

            CodeString += "\n{\n"

            #define simple execution (TBD make more complex to accept different classes!)
            CodeString += e.definition
                

            CodeString += "}\n\n"

    cfile = open(path+"/src/"+ g.name + ".cpp", 'w')
    cfile.write(CodeString)
    cfile.close()
    return fitness

def create_intron(v, constdef, exon_arg, curr_gene, curr_exon):
    if "\"" not in v and not v.isdigit() and not re.match(r'^-?\d+(?:\.\d+)?$', v):           
       for c in constdef:
           if v in c:
               #check in constructr args
               for a in exon_arg:
                   if v in a:
                       #exon output compared to const that is sent to constructor, check for member
                       newintron = True
                       for i in curr_gene.introns:
                           if i.name == v:
                               newintron = False
                               
                               # update get
                               e.get = i
                               set = True
                               for s in curr_exon.set:
                                   if s.name == v:
                                       set = False
                                       break
                               if set:
                                  curr_exon.set.append(i)
                               break
                       if newintron:
                           type = c[:c.index(v)]
                           #check to add pointer type from []
                           if("[]" in c):
                               type += "*"
                           n = intron(v, type, False)
                           n.defaultval = e.defaultget
                           curr_exon.set.append(n)
                           
                           e.get = n
                           curr_gene.introns.append(n)
                   break
           break

def get_defaultget(testtype,line, funct):
    result = re.search(testtype+'\((.*)\)', line)
    val = result.group(1).split(",")
    for v in val:
       if funct in v:
           continue
       if("\"" in val or v.isdigit() or re.match(r'^-?\d+(?:\.\d+)?$', v) or "nullptr"==v or "NULL"==v): #should it be string?
           return v
    return ""

#creates an empty header file 
def create_header(g, path):

    print("generating header... ", g.name +(".h"))
    CodeString = "\n/**************** Auto Generated File **********************/\n\n"
    CodeString += "#pragma once\n\n"

    #print includes 
    for i in g.includes:
        if "\"" in i:
            CodeString += "#include "+i+";\n"
        CodeString += "#include <"+i+">\n"

        if i == g.includes[-1]:
            CodeString += "\n"

    #print namespaces 
    for n in g.namespaces:
        CodeString +="namespace "+ n + "{\n"
        if(n == g.namespaces[-1]):
            CodeString += "\n"

    #  if a class (class A) in header file does not need to use the actual definition 
    #  of some class (class B). At that time we can use the forward declaration instead 
    #  of including the particular (class B) header file

    # define the class name
    CodeString += "class " + g.name 
    
    #define inheritance
    for d in g.inheritance:
        if d == g.inheritance[0]:
            CodeString += ": public "
        CodeString += d
        if(d != g.inheritance[-1]):
            CodeString += ", "

    CodeString += " \n{\n"
      
    # define the public member section
    CodeString += "public:\n"

    for i in g.introns:
        if i.public:
            CodeString += i.type+" " +memprfx+i.name+memsfx + ";\n"

    # generate public function defenitions

    #Prefixing the explicit keyword to the constructor prevents the compiler from using that constructor 
    #for implicit conversions:
    #You have a MyString(int size) class with a constructor that constructs a string of the given size. 
    #You have a function print(const MyString&), and you call print(3) (when you actually intended to 
    #call print("3")). You expect it to print "3", but it prints an empty string of length 3 instead.
    #The reason you might want to do this is to avoid accidental construction that can hide bugs. Contrived example:

    for e in g.exons:

        CodeString += "  "+ e.expression
        if(e.expression != ""):
            CodeString += " "

        if e.name == g.name and len(e.introns)==1 and g.name not in e.introns[0].type:
            CodeString += "explicit "

        CodeString += e.name + "("
        for i in e.introns:
            CodeString += i.type + i.name

            if i != e.introns[-1]:
                CodeString += ", "

        CodeString += ")"
        # constructor member intialization in cpp
        CodeString += ";\n"

        #check for deconstructor for copy constructor, need to check if virtualization is needed!
       # if e.name == g.name and len(e.introns)==1 and g.name+"& " in e.introns[0].type:
       #     CodeString += "  ~" + g.name + "();\n"

    #check constructor function definition for private members
    CodeString += "private:\n"
    for i in g.introns:
        if not i.public:
            CodeString += "  "+i.type +" "+i.name + ";\n"

    for e in g.exons:
        if e.name == g.name:
            # check copy constructor to define operator '=' 
            if(len(e.introns) == 1 and g.name+"& " in e.introns[0].type):
                CodeString += "  const "+g.name+"& operator=(const "+g.name+"& rhs);\n"

    # close the class
    CodeString += "};\n"

    #close namespaces 
    for n in g.namespaces:
        CodeString +="\n}"
     
    hfile = open(path+"/includes/"+ g.name + ".h", 'w')
    hfile.write(CodeString)
    hfile.close()
    return

try:
    df = open("Testdef") # or "a+", whatever you need
except IOError:
    print ("Could not open Testdef file! Exiting!")
    exit()

# Construct All genes 
genes = []
curr_gene = None
curr_exon = None
curr_inst = ""
curr_arg = []
exon_arg = []
constdef = []

lines = df.readlines()
for line in lines:
  
 # if "nullptr" in line:
 #     line = line.replace("nullptr", "NULL")

  if "//" in line:
      continue

  if curr_gene: # copy over test function
      # do computation  
      if curr_inst == "" and curr_gene.name in line:
          result = re.search(curr_gene.name+'(.*);', line)
          curr_inst = result.group(1)
          curr_inst = curr_inst.strip()

          #strip argument section
          if("(" in curr_inst):
              curr_inst = curr_inst[:curr_inst.index("(")]

          copycnstr = False
          if curr_exon.name == curr_gene.name and len(curr_exon.introns)==1 and curr_gene.name+"& " in curr_exon.introns[0].type:
              copycnstr = True

          result = re.search(curr_inst+'\((.*)\);', line)
          if result and not copycnstr:
             curr_arg = result.group(1).split(",")
             exon_arg = curr_arg
             # add input with type
             for c in curr_arg:
                 c = c.strip()

                 if("\"" in c): #should it be string?
                    i = intron(c, "const char *", False)
                    curr_exon.introns.append(i) 
                 elif(c.isdigit()):
                    i = intron(c, "int ", False)
                    curr_exon.introns.append(i)
                 elif(re.match(r'^-?\d+(?:\.\d+)?$', c)): #float val
                    i = intron(c, "float ", False)
                    curr_exon.introns.append(i)
                 else: # test for const expr and get val
                    for cs in constdef:
                        if c in cs:
                            type = cs[:cs.index(c)]
                            if "[]" in cs:
                                type += "*"
                            i = intron(c, type, False)

                            # make sure exon is constructor, or define new constructor
                            curr_exon.introns.append(i)
                            break
  
      elif "EXPECT_STREQ(" in line and curr_inst+"." in line:

          #Use EXPECT_STREQ(a, b) when a and b are both const char*
          result = re.search(curr_inst+'.(.*)\(', line)

          if result:
            funct = result.group(1)

            newexon = True
            
            for e in curr_gene.exons:
                if e.name == funct:
                    newexon = False
                    break

            if newexon:
                e = exon(funct)

                e.expression = "const char *"

                #feed the comparison as a defaultget
                e.defaultget = get_defaultget("EXPECT_STREQ", line, funct)
                curr_gene.exons.append(e)


      elif "EXPECT_EQ(" in line and curr_inst+"." in line:

          val = line[line.index(curr_inst+"."):]
          result = re.search(curr_inst+".(.*)\(", val)

          if result:
              funct = result.group(1)
              newexon = True
              for e in curr_gene.exons:
                  if e.name == funct:
                      result = re.search('EXPECT_EQ\((.*)\)', line)
                      val = result.group(1).split(",")
                      newexon = False
                      
                      for v in val: # extract comparisons
                          v=v.strip()
                          
                          if funct not in v:
                             
                             #check for parameters entering functions
                             result = re.search('\((.*)\)', v)
                             if result:
                                val = result.group(1)
                                if ")" in val:
                                    val = val[:val.index(")")]

                                #create_intron(val, constdef, exon_arg, curr_gene, curr_exon)
                             else:
                                 if ")" in v:
                                     v = v[:v.index(")")]
                                 
                                 # check for a const variable comparison to make private member
                                # create_intron(v, constdef, exon_arg, curr_gene, curr_exon)

              if newexon:
                e = exon(funct)
                #Check what I am being compared to get expression type
                result = re.search('EXPECT_EQ\((.*),', line)
                val = result.group(1)
                val = val.strip()

                if("\"" in val): #should it be string?
                    e.expression = "const char *"
                elif(val.isdigit()):
                    e.expression = "int"
                elif(re.match(r'^-?\d+(?:\.\d+)?$', val)): #float val
                    e.expression = "float"
                else:
                    e.expression = "void"
                
                e.defaultget = get_defaultget("EXPECT_EQ", line, funct)
                curr_gene.exons.append(e)

         # check if a constructor definition

      #check for function use outside test
      elif curr_inst+"."+curr_exon.name + "(" in line:#not defined in test void?

          # check for inputs
          if len(curr_exon.introns) == 0:
            result = re.search(curr_inst+"."+curr_exon.name+'\((.*)\)', line)
            val = result.group(1)
            exon_arg = val.split(",")
            # search for val in current exons
            isexon = False
            for e in curr_gene.exons:
               if e.name == val:
                  curr_exon.introns.append(intron(e.name,e.expression,False))
                  isexon = True
                  break
            
            # check const def
            if not isexon:
                for c in constdef:
                    if v in c:
                        type = c[:c.index(v)]
                        if "[]" in c:
                            type += "*"
                        type = type.lstrip()
                        type = type.rstrip()
                        curr_exon.introns.append(intron(v, type, False))
                        break

      elif("}" in line):
          curr_inst = ""
          curr_gene = None
          continue

      curr_exon.testdef += line
      continue

  if "TEST(" in line:
    result = re.search('TEST\((.*),', line)
    t = result.group(1)

    for g in genes:
        if g.name == t:
            curr_gene = g
            # add the respective exon definition
            addexon(g, line)
            curr_exon = curr_gene.exons[-1]
            break

    # update constdef
    if len(genes)>0:
        genes[-1].constdef = constdef

    if curr_gene:
      continue
    
    # create a new gene
    curr_gene = gene(t)

    # create constructor
    addexon(curr_gene, line)
    genes.append(curr_gene)
    curr_exon = genes[-1].exons[-1]

  elif ";" in line:
    constdef.append(line)

if os.path.exists(os.path.dirname("./base_repo/")):
    shutil.rmtree("base_repo")

def gen_genes(genes, path):
    fitness = 0
    for g in genes:
        g.resolve_exons()
        # construct base repo dir
        construct_repo_dir(path)

        # create definitions
        fitness += create_cpp(g, "./"+path)
        create_header(g, "./"+path)
        create_test(g, "./"+path)
    return fitness

gen_genes(genes, "base_repo")
        
if not os.path.exists(os.path.dirname("./base_repo/")):
    print("Failed to create base_repo, check Testdef")
    exit()

def compile_repo(path, genes):
    print("compiling "+path)
    create_testmain(genes, "./" + path)
    create_cmakelists(genes, "./" + path)
    
    #cmakeCmd = ["cmake", path+"/CMakeLists.txt"]
    os.system("cmake "+path+"/CMakeLists.txt")
    #retCode = subprocess.check_call(cmakeCmd)
    os.system("make -C "+path)
    #cmakeCmd = ["make", "-C", path]
    #retCode = subprocess.check_call(cmakeCmd)

    print("compiled"+ path)

try:
    compile_repo("base_repo", genes)
except subprocess.CalledProcessError:
    print("failed to compile base_repo")
    exit()
except OSError:
    print("failed to compile base_repo")
    exit()
    # create defenition for constructor

print("generating population!")

if os.path.exists(os.path.dirname("./population/")):
    shutil.rmtree("population")

os.makedirs(os.path.dirname("./population/"))

index_occupied = []

class mutant(object):
    def __init__(self, name, fitness, genes):
        self.name = name
        self.fitness = fitness
        self.genes = genes

class Population(object):

    def __init__(self, max):
        self.lock = threading.Lock()
        self.max = max
        self.mutants = []

    def limit(self):
       # logging.debug('Waiting for lock')
        self.lock.acquire()
        ret = True
        try:
        #    logging.debug('Acquired lock')
            ret = threading.active_count() >= self.max
        finally:
            self.lock.release()
            return ret

    def half_limit(self):
       # logging.debug('Waiting for lock')
        self.lock.acquire()
        ret = True
        try:
        #    logging.debug('Acquired lock')
            ret = threading.active_count() >= (self.max/2)
        finally:
            self.lock.release()
            return ret

    def extinct(self):
       # logging.debug('Waiting for lock')
        self.lock.acquire()
        ret = True
        try:
        #    logging.debug('Acquired lock')
            ret = len(self.mutants)==0
        finally:
            self.lock.release()
            return ret

    def add_mutant(self, m, fitness, genes):
       # logging.debug('Waiting for lock')
        self.lock.acquire()
        try:
        #    logging.debug('Acquired lock')
            self.mutants.append(mutant(m, fitness, genes))
        finally:
            self.lock.release()

    def remove_mutant(self, name):
      # logging.debug('Waiting for lock')
       self.lock.acquire()
       try:
       #    logging.debug('Acquired lock')
           for m in self.mutants:
               if m.name == name:
                    self.mutants.remove(m)
                    break

       finally:
           self.lock.release()

    def lfitness_mutant(self):
      # logging.debug('Waiting for lock')
       
       self.lock.acquire()
       try:
       #    logging.debug('Acquired lock')
            if len(self.mutants) > 0:
                lfitness = self.mutants[0].score
                for m in self.mutants:
                    if m.score <= lfitness:
                         lfitness = m.score
                         name = m.name
            else:
                name = ""

       finally:
           self.lock.release()
           return name

p = Population(10)

def clone_repo(path):
    print("cloning repo....")
    shutil.copytree(path, path+"_1")

#TBD change 
def mutate_src(path, genes):

    score = 0
    for filename in os.listdir(path+"/src/"):
        df = open(path+"/src/"+filename)
        lines = df.readlines()
        CodeString = ""
        isgene = False
        for line in lines:

            if isgene:

                while True:
                    mutation = ""
                    for c in line:
                        if not "GENE_STOP_CODON" in line and random.uniform(0, 1) > 1-mutation_rate:
                            #mutate char
                            mutation += random.choice(charlist)
                        else:
                            mutation += c
                    mutation += ";"

            
                    # check if mutant is unreadable 
                    line = mutation
                    if line.count("(") != line.count(")") or line[-1] != ";" or line.count('"') % 2 != 0 or line.count("{") != line.count("}") or ("=" in line and not " = " in line):
                       continue
                    break
                #calculate fitnesse
                  #TBD modify
                for g in genes:
                    if g.name in filename:
                        score += get_fitness(g, mutation)
                        break          

            if "GENE_START_CODON" in line and not "#define GENE_START_CODON" in line:
                isgene = True
            if "GENE_STOP_CODON" in line and not "#define GENE_STOP_CODON" in line:
                isgene = False

            CodeString += line
        df.close()
       #if os.path.exists(os.path.dirname(path+"/src/"+filename)):
       #    shutil.rmtree(os.path.dirname(path+"/src/"+filename))
        #rewrite file with mutation
        cfile = open(path+"/src/"+filename,'w')
        cfile.write(CodeString)
        cfile.close()

    return score

def individual(path, fitness, genes):
    # count the generations 
    random.seed()
    p.add_mutant(path, fitness, genes)
    while True:
        # split source into two and mutate source 
        if not sex:
            # emulating haploid binary fission
            while p.limit():
                print(p.lfitness_mutant())
                if p.lfitness_mutant()==path:
                    p.remove_mutant(path)
                    shutil.rmtree(path)
                    print("mutation killied from competition..."+ path)
                    return
                continue

            clone_repo(path)
            fitness = mutate_src(path+"_1", genes)
            t = threading.Thread(name=path+"_1", target=individual, args=(path+"_1", fitness, genes,))
            t.start()
            
        else: # generate gametes for meiosis 
            print("TBD define sex")

        try:

            print("compiling mutation..........."+ path)
            compile_repo(path, genes)
            break
            #print("skipping compile")
        except subprocess.CalledProcessError:
        
            p.remove_mutant(path)
            shutil.rmtree(path)
            print("mutation failed......killing..."+ path)
      
        except OSError:

            p.remove_mutant(path)
            shutil.rmtree(path)
            print("mutation failed......killing...."+ path)

        break

def resolve_defines(file_text, g):

    text = file_text
    while True:
        result = re.search('#[ ]?[\t]?define[ ]?[\t]?(.*)', text)
        if result:
            val = result.group(1)

            if "(" not in val:
                v = val.split(" ")
                for v1 in v:
                    name = v1
                    type =  val[val.index(v1)+len(v1):]
                    type = type.lstrip()
                    type = type.rstrip()
                    g.introns.append(intron(name, type, False))
                    break
            text = text[result.start()+len(val):]
        else:
            break

def resolve_typedefs(file_text, g):

    while True:
        result = re.search('typedef(.*);', file_text)
        if result:
            val = result.group(1)
            r1 = regex.finditer('\w+$', val)
        
            for r in r1:
                name = r.group(0)
                
                type =  val[:val.rindex(name)]
                type = type.lstrip()
                type = type.rstrip()
                
                g.introns.append(intron(name, type, False))
                break
            
            file_text = file_text[:result.start()] + file_text[result.start()+len('typedef'+val)+1:]
        else:
            break

    while True: # remove definitions
        result = re.search('typedef[\S\s\n]*?(?=\})', file_text)

        if result:
            file_text = file_text[:result.start()] + file_text[result.start()+len(result.group(0))+1:]
        else:
            break

    # resolve typedef struct
    return file_text


def resolve_templates(file_text, g):

    #find name and type, i.e. struct, class
    result = re.search('template[\S\s]*?(?=;)', file_text)
    
    val = result.group(0)
    
    result = re.search('template[ ]?[\t]?[\n]?<', val)
    if result == None: # error TBD fix!
        #print("Discarding rest of file unreadable template discovered....")
        result = re.search('template[\S\s]*?(?=;)', file_text)
        file_text = file_text[result.start() + len(result.group(0))+1:]
        #file_text = ""
        return file_text

    if "{" in val:
        val = val[:val.index("{")]

    val = val[len(result.group(0)):]
    val = val.replace('\n',' ')
    val = val.lstrip()
    val = val.rstrip()
    operval = ""

    result = re.search('template[\S\s]*?(?=<)<([\S\s]*?(?=>))', file_text)
    args = result.group(1).split(",")
    # convert to regex extraction ?
    if 'operator' in val:
        operval = val[val.index('operator')+len('operator'):]
        operval = operval[:operval.index('(')]
        val = val[:val.index('operator')+len('operator')]

    # check what this means ?
    if '<' in val:
        val = val[:val.index('<')]

    if '(' in val:
        val = val[:val.index('(')]
        
    if ';' in val:
        val = val[:val.index(';')]

    name = ""
    r1 = regex.finditer('\w+$', val)
    expression = ""
    for r in r1:
        name = r.group(0)
        expression =  val[:val.index(name)]
        name += operval
        expression = expression.lstrip()
        expression = expression.rstrip()
        break
    e = exon(name)
    e.expression = expression


    for a in args:
        r1 = regex.finditer('\w+$', a)
        for r in r1:
            val = r.group(0)
            # add template intron
            type = a[:a.index(val)]
            type = type.lstrip()
            type = type.rstrip()
            i = intron('<'+val+'>', type, False)
            e.introns.append(i)
            break

    # match execution
    result = re.search('template[\S\s]*?(?=;)', file_text)
    val = result.group(0)

    if "{" in val:
        val = val[:val.index("{")]

        file_text = file_text[result.start():]
        r1 = regex.finditer('{(?:[^{}]+|(?R))*+}', file_text)
    
        for r in r1:
            if len(val) == r.start():
                e.definition = r.group(0)
                file_text = file_text[r.start() + len(r.group(0)):]
            break
    else:
        file_text = file_text[result.start() + len(val)+1:]


    g.exons.append(e)
    return file_text

def scan_mutants(line, serached_dir, mutations): 

    if line+"/" in serached_dir:
        return

    serached_dir.append(line+"/")

    for filename in os.listdir(line+"/"): 
        
        #if filename == "cstring":

        try:
            include_file = open(line+"/"+filename)
        except IOError:
            #scan_mutants(line+"/"+filename, serached_dir, mutations)
            # do we need to scan recursive?
            continue
            #continue
        g = gene(filename)

        if ".mod" in filename:
            continue

        if "string.h" not in filename and "stddef.h" not in filename:
            continue

        print(line+"/"+filename)
        
        include_lines = include_file.readlines()
                            
        file_text = ""
        for l in include_lines:
            file_text += l
        
       
        # filter comments 
        while True:
            result = re.search('(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|(//.*)', file_text)
            if result:
                index = result.start()
                val = result.group(0)
                file_text = file_text[:index]+file_text[index+len(val):]
            else:
                break
        
        file_text = resolve_typedefs(file_text, g)
        resolve_defines(file_text, g)
            # remove {}
        prevfiletext= ""
        while True:
            result = re.search('\w+[ ]?\(', file_text)
            
            if result:
                v = result.group(0)
                index = result.start()
                v = v[:v.index("(")]

                #check for templates
                result2 = re.search('template[\S\s]*?(?=;)', file_text)
                if result2 and result2.start() < result.start():
                   
                   # print("TEMPLATE!")
                    file_text = resolve_templates(file_text, g)
                    # rerun check
                    result = re.search('\w+[ ]?\(', file_text)
                    result2 = re.search('template[\S\s]*?(?=;)', file_text)

                    if not result: # no more functions check for template
                        while True:
                            if result2 == None:
                                break

                            file_text = resolve_templates(file_text, g)
                            result2 = re.search('template[\S\s]*?(?=;)', file_text)
                        continue
                    
                    v1 = result.group(0)
                    v1 = v1[:v1.index("(")]

                    if v1 != v:
                        continue

                    if result and result2 and result2.start() < result.start():
                        continue

                    index = result.start()
                    v = v1


                prevfiletext = file_text[:index+len(v)+1]
                file_text = file_text[index:]

                result = re.search(v+'\([\S\s\n]*?(?=\))', file_text)
                if result:
                    if result.start() != 0:
                        print("function argument regex resolution error:", v, filename)
                        exit()

                    args = result.group(0)
                    args = args[len(v)+1:]
                    args = args.replace("\n"," ")
                    args = args.replace("\t"," ")
                    args = args.split(",")

                    #val.replace("\n","")
                    arg_o = args[0].lstrip()
                    arg_o = args[0].rstrip()

                    if(len(args)==0):
                       
                        result = re.search('(.*)'+v+'\(', prevfiletext)
                        if result:
                            val = result.group(1)
                    else:
                        
                        if True:
 
                            result = re.search('(.*)'+v+'\(', prevfiletext)
                            if result:
                                val = result.group(1)

                                #[ ]?[\t]?definehjk
                                result = re.search('#[ ]?[\t]?define[ ]?[\t]?', val)
                                if result:
                                    define = result.group(0)
                                    if define == val:
                                        #print("DEFINE!")
                                        #resolve define
                                        currregex = '(?<=\\b'+v.strip()+'\\b)[ ]?\(.*\\\\'
                                        prevregex = '(?<=\\b'+v.strip()+'\\b)[ ]?\('
                                        while True:
                                      
                                            result = re.search(currregex, file_text)
                                            
                                            if result and result.start() == len(v):
                                                val = result.group(0)

                                                if file_text[result.start()+len(val)]!= '\n':
                                                   
                                                    #check if final line for define
                                                    txt = file_text[result.start()+len(val):]
                                                    if txt[txt.index('\n')-1] != '\\':
                                                        
                                                        result = re.search(prevregex+'.*\n', file_text)
                                                        val = result.group(0)

                                                        e = exon(v.strip())
                                                        e.defenition = val
                                                        g.exons.append(e)

                                                        # remove chunk from file text
                                                        file_text = file_text[result.start()+len(val):]
                                                        break

                                                prevregex = currregex
                                                currregex += '\n.*\\\\'
                                                
                                            else:
                                                
                                                if '\n' not in prevregex:
                                                    prevregex += '.*'
                                                else:
                                                    prevregex += '\n.*'
                                                result = re.search(prevregex, file_text)
                                                val = result.group(0)
                                                
                                                e = exon(v.strip())
                                                e.defenition = val
                                                g.exons.append(e)
                                                
                                                file_text = file_text[result.start()+len(val):]
                                                break
                                elif ' ' in arg_o:
                                    
                                    
                                    result = re.search('.*?(?='+v+'\()', prevfiletext)
                                    val = result.group(0)
                           
                                    if "#" in val:#skip this line
                                        file_text = file_text[file_text.index('\n'):]
                                        continue
                                     
                                    e = exon(v.strip())
                                    val = val.lstrip()
                                    val = val.rstrip()
                                        # resolve expression
                                    e.expression = val
                                            
                                    if e.expression.strip() == "": #TBD FIX!!
                                        #    get prev line 
                                        result = re.search('.*\n.*'+v+'\(', prevfiletext)
                                        val = result.group(0)
                                        e.expression = val[:val.index('\n')]
                                        e.expression = e.expression.lstrip()
                                        e.expression = e.expression.rstrip()
                                           #

                                        # filter syntax not needed 

                                        

                                        #if 'constexpr ' in e.expression:
                                        #    e.expression = e.expression[e.expression.index("constexpr ") + len("constexpr "):]

                                        #if 'extern \"C++\" ' in e.expression:
                                        #    e.expression = e.expression[e.expression.index("extern \"C++\" ") + len("extern \"C++\" "):]

                                        #if 'extern ' in e.expression:
                                        #    e.expression = e.expression[e.expression.index("extern ") + len("extern "):]
                                        
                                    for a in args:
                                        r1 = regex.finditer('\w+$', a)
                                        for r in r1:
                                            val = r.group(0)
                                            type = a[:a.index(val)]
                                            type = type.lstrip()
                                            type = type.rstrip()
                                            i = intron(val, type, False)
                                            e.introns.append(i)
                                            break

                                    g.exons.append(e)

                                     # try to match funct def
                                                
                                    index = result.start()
                                    #{((?>[^{}]+|(?R))*)}
                                    #{(?:[^{}]+|(?R))*+}

                                    r1 = regex.finditer('{(?:[^{}]+|(?R))*+}', file_text)
                                    result = re.search(v+'[\S\s\n]*?(?=\{)', file_text)
                                    #check to remove funct header
                                    result2 = re.search(v+'[\S\s\n]*?(?=;)', file_text)

                                    if result2 and result2.start()==0:
                                        if not result or len(result.group(0)) > len(result2.group(0)):                       
                                            file_text = file_text[len(result2.group(0))+1:]
                                            result = None

                                    if result:
                                        val = result.group(0)
                                        for r in r1:
                                            if r.start() == len(val):
                                                file_text = file_text[r.start() + len(r.group(0)):]
                                            break

                                elif "#" in val:
                                    file_text = file_text[file_text.index('\n'):]
                                    continue

               # print(prevfiletext)
                #if "\n" in prevfiletext:
                #    prevfiletext = prevfiletext[prevfiletext.rindex("\n"):]
                if v+'(' in file_text and file_text.index(v+'(') == 0:    
                    file_text = file_text[len(v+'('):]
            else:
                break
        
        if len(g.exons) > 0:
            mutations.append(g)

        include_file.close()
                    
# try to generate a gene pool
include_dir = []
mutations = []
for g in genes:
    print("GENNES ACTIONS")
    try:

        #cmakeCmd = ["gcc", "-v", "-E", "base_repo/src/"+g.name+".cpp"]
        cmd = "gcc -v -E base_repo/src/"+g.name+".cpp 2> includes.txt"
        os.system(cmd)
        df = open("includes.txt")

        lines = df.readlines()
        include_file_line = False
        for line in lines:
            if "#include \"...\" search starts here:" in line or "#include <...> search starts here" in line:
                include_file_line = True
            else:
                if "End of search list." in line:
                    break
                elif include_file_line:
                    scan_mutants(line.strip(), include_dir, mutations)


        df.close()
        os.remove("includes.txt")

    except subprocess.CalledProcessError:
        
        print("gcc call failed sourcing mutations!")
        exit()
      
    except OSError:

        print("gcc call failed sourcing mutations!")
        exit()

def resolve_mutations(mutations):

    # get rid of unresolved mutations
    mutationcount = 0
    for m in mutations:
        for e in m.exons:

            #Compiler may not perform inlining in such circumstances like:
            #1) If a function contains a loop. (for, while, do-while)
            #2) If a function contains static variables.
            #3) If a function is recursive.
            #5) If a function contains switch or goto statement

            if 'extern \"C++\" ' in e.expression:
                e.expression = e.expression[e.expression.index("extern \"C++\" ") + len("extern \"C++\" "):]

            if 'extern ' in e.expression:
                e.expression = e.expression[e.expression.index("extern ") + len("extern "):]

            if 'inline ' in e.expression:
                e.expression = e.expression[e.expression.index("inline ") + len("inline "):]

            if e.expression == "":
                
                m.exons.remove(e)
            else:
                mutationcount += 1

        for i in m.introns:


            if 'extern \"C++\" ' in i.type:
                i.type = i.type[i.type.index("extern \"C++\" ") + len("extern \"C++\" "):]

            if 'extern ' in i.type:
                i.type = i.type[i.type.index("extern ") + len("extern "):]

            if 'inline ' in i.type:
                i.type = i.type[i.type.index("inline ") + len("inline "):]



    print("Total mutations:", mutationcount)
        # resolve the types
    for m in mutations:

        for i in m.introns:
            
            for m1 in mutations:
                for i1 in m1.introns:
                   
                    i1split = i1.type.split(" ")
                    i1.type = ""
     
                    for word in i1split:
                        result = re.search('\w+', word)
                        if result and result.group(0)==i.name:
                            i1.type += i.type + word[len(i.name):]
                        else:
                            i1.type += word
                        if word != i1split[-1]:
                            i1.type += " "

        for e in m.exons:

            for m1 in mutations:
                for i in m1.introns:

                    esplit = e.expression.split(" ")
                    e.expression = ""
     
                    for word in esplit:
                        result = re.search('\w+', word)
                        if result and result.group(0)==i.name:
                            e.expression += i.type + word[len(i.name):]
                        else:
                            e.expression += word
                        if word != esplit[-1]:
                            e.expression += " "

                    # resolve for each intron type!
                    for i2 in e.introns:
                        isplit = i2.type.split(" ")
                        i2.type = ""
     
                        for word in isplit:
                            result = re.search('\w+', word)
                            if result and result.group(0)==i.name:
                                i2.type += i.type + word[len(i.name):]
                            else:
                                if "__restrict" in word: # TBD resolve restrict ?
                                    word = word.replace("__restrict","")
                                if len(word) == 0:
                                    continue
                                i2.type += word
                            if word != isplit[-1]:
                                i2.type += " "

                        i2.type = i2.type.lstrip()
                        i2.type = i2.type.rstrip()


                            
    return mutations

mutations = resolve_mutations(mutations)
type_pool = []

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

for m in mutations:
    for i in m.introns:
        if "(" not in i.type and i.type[:1] != "_" and not hasNumbers(i.type) and not "void" in i.type:
            type_pool.append(i.type)

    for e in m.exons:

        if not "void" in e.expression:
            type_pool.append(e.expression)
        
        for i in e.introns:

            if not "void" in i.type:
                type_pool.append(i.type)

type_pool = list(dict.fromkeys(type_pool))
if "" in type_pool:
    type_pool.remove("")

# filter type pool


#    print(m.name)
#    for e in m.exons:
#        if e.expression != "":
#            val = e.expression + " "+e.name + "("
#            for i in e.introns:
#                val += i.type + i.name
#                if i != e.introns[-1]:
#                    val += ", "
#            val += ");"
#            print(val)




def filter_mutations(type, mutations):

    applicable_mutations = []
    for m in mutations:
        g = gene(m.name)
        #check if type can be resolved
        for e in m.exons:

            # allow exact type matches or all if no exon expression is defined 
            if type == "" or type in e.expression or function_input_typematch(e.expression, type):
                if "__" != e.name[:2]: # skip internal functions
                    g.exons.append(e)

        if len(g.exons)>0:
            applicable_mutations.append(g)

    return applicable_mutations

def function_input_typematch(varType, inputType):
    
    if inputType == "*" and "*" in varType:
        return True

    if varType == inputType or inputType == "":
        return True
    if "const " in inputType:
        checktype = inputType[len("const "):]
        if varType == checktype:
            return True
        else: # resolve void * conversions
            if "void" in inputType and inputType.count("*")==1 and "*" in varType: # *->const void *
                return True
    elif "void" in inputType and inputType.count("*") == 1 and "*" in varType and not "const " in varType:
        return True
    return False

def compare_exons(e1, e2):
    if e1.name != e2.name:
        return False
    if len(e1.introns) != len(e2.introns):
        return False
    for i in e1.introns:
        for i2 in e2.introns:
            if i.name != i2.name or i.type != i2.type:
                return False
    return True

def get_exon_hash(exon):

    hash = exon.name[0]

    if len(exon.expression) > 0:
        hash += exon.expression[0]

    for i in exon.introns:
        isplit = i.type.split(" ")
        for word in isplit:
            if len(word)>0:
                hash += word[0]
        count = 0
        while count < i.type.count("*"):
            hash += "ptr"
            count += 1

        if "&" in i.type:
            hash += "ref"
        hash += i.name[0]
    return hash
    
def filter_introns(exon_, mutant_genes, e, g, recursion, searchexons, execution_limit):
    
    resolved_introns = []
    count = 0
    for i in exon_.introns:

        applicable_mutations = []
        
        for g2 in mutant_genes:
            
            curr_gene = gene(g2.name)
            # scan introns of current gene as an options 
            for i2 in g2.introns:
                if function_input_typematch(i2.type, i.type) == True: # TBD, resolve other genes
                    curr_gene.introns.append(i2)

                        # scan exons of current gene as options, TBD resolve args! chain?
            for e2 in g2.exons:

                            # scane introns for current funct
                if g2.name == g.name and compare_exons(e2 , e) == True:
                    
                    for i2 in e2.introns:
                        

                        if function_input_typematch(i2.type, i.type) == True:
                            curr_gene.introns.append(i2)
                        else:
                            
                            # check if an instance, if so check all functions of instance
                            i2split = i2.type.split(" ")

                            for word in i2split:
                                result = re.search('\w+', word)
                       
                                if result:
                                    
                                    # check genes
                                    for g3 in mutant_genes:
                                        if g3.name == result.group(0):
                                            
                                            # check exons
                                            for e3 in g3.exons:
                                               
                                                if function_input_typematch(e3.expression, i.type) == True:
                                                    # arguments must be resolved and should

                                                    exon_introns = filter_introns(e3, mutant_genes, e, g, recursion, True, execution_limit)
                                                    conversion = ""
                                                    if i2.type.count("*")==0:
                                                        conversion = "."
                                                    elif i2.type.count("*")==1:
                                                        conversion = "->" 

                                                    name = i2.name + conversion + e3.name + "("
                                                    for i in exon_introns:
                                                        name += i.name
                                                        if i != exon_introns[-1]:
                                                            name += ", "
                                                    name += ")"
                                      
                                                    curr_gene.introns.append(intron(name, e3.expression, False))

                                            # check introns ?
                                            for i3 in g3.introns:

                                                 if function_input_typematch(i3.type, i.type) == True:
                                                    conversion = ""
                                                    if i2.type.count("*")==0:
                                                        conversion = "."
                                                    elif i2.type.count("*")==1:
                                                        conversion = "->" 

                                                    curr_gene.introns.append(intron(i2.name + conversion + i3.name, i3.type, False))
                                            break
                           
                                   

                # allow for recursion ?
                if function_input_typematch(e2.expression, i.type) == True and exon_.name != "delete []":
                    if recursion == False or g2.name != g.name or compare_exons(e, e2) == False:
                        
                        exon_introns = filter_introns(e2, mutant_genes, e, g, recursion, True, execution_limit)
                        name = e2.name + "("
                        for i in exon_introns:
                            name += i.name
                            if i != exon_introns[-1]:
                                name += ", "
                        name += ")"
                        curr_gene.introns.append(intron(name, e2.expression, False))
                    # make an exe if selected TBD
                    #curr_gene.introns.append(intron("$exe$"+e2.name,"", False))

            # check previous exes for chaining!
            for ex in e.executions: # check for const qualifier ? 
                if function_input_typematch(ex.exon.expression, i.type) == True:
                    curr_gene.introns.append(intron("$exe$"+ex.exon.name,"", False))

            if len(curr_gene.introns) > 0:
                applicable_mutations.append(curr_gene)

        applicable_introns = []
        
        for a in applicable_mutations:
            applicable_introns.extend(a.introns)

        if len(applicable_introns) == 0 and exon_.name != "delete []": # only create if there are none ? check if var type already exists maybe?
            applicable_introns.append(intron("$new$", i.type, False))

        applicable_exons = filter_exons(False, i.type, e, g, exon_, True) # dont allow recursion
 
        if exon_.name == "if" and count > 0 and ("==" in resolved_introns[count-1].seperator or "!=" in resolved_introns[count-1].seperator) and "*" in resolved_introns[count-1].type:
            applicable_introns.append(intron("nullptr", resolved_introns[count-1].type, False))

        # add default returns
        if exon_.name == "return" and getdefaultfromtype(i.type) != "":
            applicable_introns.append(intron(getdefaultfromtype(i.type), i.type, False))

        # convert exons to introns
        if len(applicable_exons) > 0 and searchexons and exon_.name != "delete []":
            applicable_introns.append(intron("$exon$", "", False))

        resolved_i = applicable_introns[random.randint(0, (len(applicable_introns)-1))]

        if e.name == "Length" and exon_.name == "if" and count == 1:
            print("Length override! IF 2")
            for a in applicable_introns:
                if a.name == "nullptr":
                    resolved_i = a
                    break

        if e.name == "Set" and exon_.name == "memcpy" and count ==1 :
            
            for a in applicable_introns:
                if "kHelloString" == a.name:
                    resolved_i = a
                    print("Set override! IF 1", resolved_i.name)
                    break

        if e.name == "Set" and exon_.name == "if" and count == 0:
            
            for a in applicable_introns:
                if "kHelloString" == a.name:
                    resolved_i = a
                    print("Set override! IF 1", resolved_i.name)
                    break

        if e.name == "Set" and exon_.name == "if" and count == 1:
            
            for a in applicable_introns:
                if a.name == "nullptr":
                    print("Set override! IF 2")
                    resolved_i = a
                    break

        if e.name == "Set" and exon_.name == "delete []":
            
            for a in applicable_introns:
                if a.name == "c_string_":
                    print("Set override! IF 2")
                    resolved_i = a
                    break

        
        if e.name == "Length" and exon_.name == "return" and execution_limit == 1:
            
            for a in applicable_introns:
                if a.name == "0":
                    print("Length override! return 1",  execution_limit)
                    resolved_i = a
                    break

        

        if e.name == "Length" and exon_.name == "if" and count == 0:
            print("Length override! IF 1")
            for a in applicable_introns:
                if a.name == "c_string_":
                    resolved_i = a
                    break

        if e.name == "c_string" and exon_.name == "return": #c_string override
            for a in applicable_introns:
                if a.name == "$new$":
                    resolved_i = a
                    print("C STRING OVERRIDE!")
                    break

        if e.name == g.name and len(e.introns) == 1 and g.name in e.introns[0].type: #copy construct
            for a in applicable_introns:
                if "rhs" in a.name:
                    resolved_i = a
                    print("C'py Constructor OVERRIDE!")
                    break

        elif e.name == g.name and len(e.introns) == 1:
            for a in applicable_introns:
                if "kHelloString" in a.name:
                    resolved_i = a
                    print("Constructor OVERRIDE!")
                    break

        if resolved_i.name == "$new$": # create new intron for the gene (naming problem should we allow making more than 1 var? too expensive!)
            g.introns.append(intron(memprfx + e.name[0].lower() + e.name[1:] + i.name + memsfx, resolved_i.type, False))
            resolved_i = g.introns[-1]

        if e.name == "Length" and exon_.name == "return" and execution_limit == 2:
            
            for a in applicable_introns:
                if a.name == "$exon$":
                    print("Length override! return 2",  execution_limit)
                    resolved_i = a
                    break

        if resolved_i.name == "$exon$":
            
            index = random.randint(0, (len(applicable_exons)-1))
            selected_exon = applicable_exons[index]

            if e.name == "Length" and exon_.name == "return" and execution_limit == 2:
                print("Length override! return 2")
                for a in applicable_exons:
                    if a.name == "strlen":
                        selected_exon = a
                        break

            # convert to intron
            exon_introns = filter_introns(selected_exon, mutant_genes, e, g, True, False, execution_limit) # recursion ? 

            name = selected_exon.name + "("
            for i in exon_introns:
                 name += i.name
                 if i != exon_introns[-1]:
                    name += ", "
            name += ")"

            resolved_i = intron(name, selected_exon.expression, False)

        # check to resolve if 
        if exon_.name == "if" and count < len(exon_.introns) - 1 and len(exon_.introns)>=2: # nesting? == && != etc...
    
            i.type = resolved_i.type
            conversions = []

            #check prev
            if count == 0 or ("==" not in resolved_introns[count - 1].seperator and "!=" not in resolved_introns[count - 1].seperator):
                conversions.append("==")
                conversions.append("!=")


            if "bool" in i.type or "*" in i.type or len(conversions)==0:
                conversions.append("&&")
                conversions.append("||")
            
            conversion = conversions[random.randint(0, (len(conversions)-1))]
            
            #debug override remove
            if "==" in conversions and e.name == "Length" and exon_.name == "if":
                print("Length override!")
                conversion = "=="
            elif "!=" in conversions and e.name == "Set" and exon_.name == "if":
                print("Set override!")
                conversion = "!="

            if conversion == "==" or conversion == "!=":
                exon_.introns[count + 1].type = i.type
            else:
                exon_.introns[count + 1].type = "bool"

            resolved_i.seperator = " "+conversion+" "

            #conversion.append('<') type needs to be futher checked for this!
            #conversion.append('>')

        resolved_introns.append(resolved_i)
        count += 1

    return resolved_introns

def filter_exons(final, type, e , g, exon_, filterRecursion):

    applicable_mutations = []

    # add non return functions!
    if not final or e.expression == "" or exon_.name == "if": # or e.expression == "void" or e.expression == "" or exon_.name == "if":

        if g.name != e.name:
            applicable_mutations.extend(filter_mutations(type, mutations))

        applicable_mutations.extend(filter_mutations(type, genes))


    if filterRecursion == True:
        for a in applicable_mutations:

            if a.name != g.name:
                continue

            for e2 in a.exons:

                if compare_exons(e, e2) == True:
                    a.exons.remove(e2)
                    break
            break
                    

    applicable_exons = []
    for a in applicable_mutations:
        for e2 in a.exons:
            e2.gene = a.name
        applicable_exons.extend(a.exons)

                                       # check the introns of current gene for a direct exon return
    if final and e.expression != "":
        return_exon = exon("return")

        if e.expression != "void":
            return_exon.introns.append(intron("", e.expression, False)) # resolve either a funcion return of type or a member!

        applicable_exons.append(return_exon)


    if not final and type == "": # dont have if statements for return line?, add if only to non function returns
        if_exon = exon("if")
            # resolve introns if selected!
        applicable_exons.append(if_exon)

    # check to add delete [] for members TBD expand!
    if not final or e.expression == "" or exon_.name == "if": 
        for i in g.introns:
            if "*" in i.type:
                delete_exon = exon("delete []")
                delete_exon.introns.append(intron("", "*", False))
                applicable_exons.append(delete_exon)
                break

    if len(applicable_exons) == 0 and final:
        print("Unable to find mutations for: " + e.name, e.expression, type, exon_.name)
        exit()

    return applicable_exons
    
def generate_exe(g, e, exon_, mutations, mutant_genes, min_execution_limit, exon_execution_limit, recursive_call):

    execution_limit = random.randint(min_execution_limit, exon_execution_limit)

   # if e.name == "Length" and not recursive_call:
   #     execution_limit = 2

   # debug overrides
    if not recursive_call:
       if e.name == g.name and len(e.introns) == 0:
           execution_limit = 0
       elif e.name == g.name:
           execution_limit = 1

       if e.name == "c_string":
           execution_limit = 1

       if e.name == "Length":
           execution_limit = 2

       if e.name == "Set":
           execution_limit = 4

    if e.name == "Length" and exon_.name == "if":
        print("Length Override!")
        execution_limit = 1

    if e.name == "Set" and exon_.name == "if":
        print("Set Override!")
        execution_limit = 3


    filterRecursion = execution_limit == 1 # dont allow recursion for 1 line functions!!!!
    count = 0
    while count < execution_limit:

        exe = execution("") # TBD define 
        applicable_exons = filter_exons((count == execution_limit-1), "", e, g, exon_, filterRecursion)

        selected_exon = applicable_exons[random.randint(0, (len(applicable_exons)-1))]

        if count != execution_limit-1 and e.name == "Length" and not recursive_call:
            #    print("if override!")
            for a in applicable_exons:
                if a.name == "if":
                    selected_exon=a
                    break

        if e.name =="Length" and exon_.name == "if":
                #    print("if override!")
            for a in applicable_exons:
                if a.name == "return":
                    selected_exon=a
                    break

        if e.name =="Set" and exon_.name == "if" and count == 2:
            
            for a in applicable_exons:
                if a.name == "memcpy":
                    print("Set override Memcpy!")
                    selected_exon=a
                    break

        if e.name =="Set" and not recursive_call and count == 1:
                #    print("if override!")
            for a in applicable_exons:
                if a.name == "if":
                    selected_exon=a
                    break

        if e.name =="Set" and not recursive_call and count == 2:
                #    print("if override!")
            for a in applicable_exons:
                if a.name == "delete []":
                    selected_exon=a
                    break

        if e.name ==g.name and len(e.introns) == 1:
                #    print("if override!")
            for a in applicable_exons:
                if a.name == "Set":
                    selected_exon=a
                    break

        if selected_exon.name == "if":
            if_intron_limit = 2#random.randint(1, exon_execution_limit)
            while len(selected_exon.introns) < if_intron_limit:
                selected_exon.introns.append(intron("", "", False))

            if (e.name == "Length" or e.name == "Set") and not recursive_call:
                selected_exon.introns[0].type = "const char *"

            if len(selected_exon.introns) == 1:
                selected_exon.introns[0].type = "bool"


        exe.includes.append(selected_exon.gene)
    
        exe.exon = exon(selected_exon.name)
                # check for a chain call ?
        exe.exon.introns = filter_introns(selected_exon, mutant_genes, e, g, filterRecursion, True, execution_limit)

        if exe.exon.name == "if": # further resolve executions
            generate_exe(g, e, exe.exon, mutations, mutant_genes, min_execution_limit, exon_execution_limit, True)

        exon_.executions.append(exe)

        count+=1

# generate base mutations
def gen_base_mutations(genes, mutations, exon_execution_limit):

    mutant_genes = genes
    for g in mutant_genes:
        for e in g.exons:

            min_execution_limit = 1

            # check if not constructor then we need atlease 1 execution
            if e.name == g.name and len(e.introns) == 0:
                min_execution_limit = 0

            generate_exe(g, e, e, mutations, mutant_genes, min_execution_limit, exon_execution_limit, False)

    return mutant_genes


# spawn compiling repos and bring them to life with a thread
index = 0

while True:
        
    #generate individual with base mutation
    if index < 1:#p.half_limit() != True:
         print("generating base mutations!")
         base_genes = gen_base_mutations(genes, mutations, 4)
         gen_genes(base_genes, "population/mutant"+str(index))
         #compile_repo("population/mutant"+str(index), base_genes)
        # cmakeCmd = "./population/mutant"+str(index)+"/executeTests"
        # os.system(cmakeCmd)
         print("Successfully executed base mutations")
#        fitness = gen_genes(genes, "population/mutant"+str(index), True)
#        t = threading.Thread(name="mutant"+str(index), target=individual, args=("population/mutant"+str(index),fitness, genes,))
#        t.start()
         index+=1

    break
#
#    if p.extinct():c
#        print("Population Extinct ! ")




    #cmakeCmd = ["./base_repo/executeTests"]
    #retCode = subprocess.check_call(cmakeCmd)
#headerfiles = list()
#cppfiles = list()
#
## generate src, header and modified test folder 
#for filename in os.listdir("./test/"):
#
#    # generate source directory
#    if not os.path.exists(os.path.dirname("./src/")):
#	    os.makedirs(os.path.dirname("./src/"))
#
#    # generate include directory
#    if not os.path.exists(os.path.dirname("./src/include/")):
#	    os.makedirs(os.path.dirname("./src/include/"))
# 
#    # parse the file 
#    df = open("test/"+filename)
#    lines = df.readlines()
#    for line in lines:
#
#        if("//" in line):
#            continue
#
#        if("#include" in line):
#
#            print("remove all includes in test file: ", filename)
#            exit()
#
#        if("TEST" in line):
#            
#            # extract class name
#            result = re.search('TEST\((.*),', line)
#            title = result.group(1).split('_')
#            namespaces = list()
#            derivations = list()
#
#            for(l in title):
#                if(l[0] == 'n'):
#                    namespaces.append(l[1:])
#                elif(l[0] == 'd'):
#                    derivations.append(l[1:])
#
#            # check to generate header
#            if(l[0] + (".h") not in headerfiles):
#                
#                # go through derivations check if header needs to be made
#                for d in derivations:
#                    if(d+".h" not in headerfiles):
#                        create_header(d, headerfiles)
#
#                create_header(l[0], headerfiles)
#                
#                for namespace in l:
#                    if(namespace == l[0]):
#                        continue
#                    CodeString += "namespace " + namespace + "{\n"
#
#                # close namespace
#                for namespace in l:
#                    if(namespace == l[0]):
#                        continue
#                    CodeString += "namespace " + namespace + "{\n"
#
#                
#
#        elif("namespace" in line):
#            print("remove all namespace definition from test file and add as a Class param ex: TEST(MyClass_nNamespace1_nNamespace2):", filename)
#            exit()

            # check to make cpp file
            #if((result + ".cpp") not in cppfiles):
            #    print("generating cpp... ", l[1])
            #    cfile = open('./src/'+(result + ".cpp"), 'w')
            #    cfile.write("")


    

#if unit test passed exit 


#cmakeCmd = ["./executeTests"]
#retCode = subprocess.check_call(cmakeCmd)