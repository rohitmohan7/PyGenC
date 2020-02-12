
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

# represents a gene function
class exon(object):

    def __init__(self, name):
        self.name = name
        self.introns = []
        self.sub_exons = [] # if present must be istantiated or used via static i.e. 'gene::' or 'gene.' 
        self.expression = ""
        self.testdef = ""
        self.get = None # defines exon direct intron get
        self.defaultget = ""
        self.set = [] # defines exon direct intron sets
        self.public = False
        self.definition = ""
        self.executions = []

class intron(object):

    def __init__(self, name, type, public):
        self.name = name
        self.type = type
        self.public = public
        self.defaultval = ""

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
            e.definition = ""
            continue
            for s in e.set:
                # check if constructor
                if e.name == self.name:
                    #check for another function that may be setting this 
                    for e2 in g.exons:
                        set = False
                        if e2.name != self.name:
                            for s2 in e2.set:
                                if s2.name == s.name:
                                    e.definition += "   "+e2.name+"("+s.name+");\n"
                                    set = True
                                    break
                            if set:
                                break
                else:
                     if "const char *" in s.type:
                        print("nothing")
                        #CodeString += "  if ("+s.name+")\n  {\n"
                        #CodeString += "    const size_t len = strlen("+s.name + ");\n"
                        #CodeString += "    char* const clone = new char[len + 1];\n"
                        #CodeString += "    memcpy(clone, "+s.name+", len + 1);\n"
                        #CodeString += "\n    if ("+memprfx + s.name + memsfx +" != NULL)\n    {\n"
                        #CodeString += "      delete[] "+memprfx + s.name + memsfx + ";\n    }\n\n"
                        #CodeString += "    "+memprfx + s.name + memsfx +" = clone;\n" 
                        #CodeString += "  }\n  else\n  {\n"
                        #CodeString += "    delete[] "+memprfx + s.name + memsfx + ";\n"
                        #CodeString += "    "+memprfx + s.name + memsfx + " = NULL;\n  }\n"

                        #add to includes 
                        #if not "string" in g.includes:
                        #    g.includes.append("string")

            if e.get:

                for i in self.introns:
                    if i.name == e.get.name:
                        # check return type to return member here within all C++ rules and make a sensible conversion
                        if e.expression in i.type:
                            e.definition += "   return "+memprfx+e.get.name+memsfx+";\n"
                        #elif mutate:
                        #    CodeString += "GENE_START_CODON\n"
                            # add random mutation 
                        #    mutation = gen_mutation(g)
                            
                          #  fitness += get_fitness(g, line)
                        #    CodeString += mutation
                            
                        #    CodeString += "\nGENE_STOP_CODON\n"

                        #elif "int" in e.expression and "const char *" in i.type: # send length
                            #CodeString += "   if(!"+memprfx+i.name+memsfx+")\n    return "+e.defaultget+";\n"
                            #CodeString += "   return strlen("+memprfx+i.name+memsfx+");\n"
                        #    if not "string" in g.includes:
                        #        g.includes.append("string")
                        break

            # define copy constructor
            if(e.name == self.name and len(e.introns) == 1 and self.name+"& " in e.introns[0].type):            
                # find a method to copy from rhs all introns
                for i in self.introns:
                    method = False
                    for e2 in self.exons:
                        if e2.name != g.name:
                            for s in e2.set:
                                if s.name == i.name:
                                    e.definition += "  "+e2.name+"("+e.introns[0].name+"."+memprfx + i.name + memsfx+");\n"
                                    method = True
                                    break
                        if method:
                            break
                    if not method:
                        e.definition += "  " +memprfx + i.name + memsfx+ " = "+e.introns[0].name +"."+memprfx + i.name + memsfx+";\n"

                e.definition += "}\n\n"

                # define deconstructor
                e.definition += self.name + "::~"+self.name+"()\n{\n"
                for i in self.introns:
                    if "*" in i.type:
                        e.definition += "  delete[] "+memprfx+i.name+memsfx+";\n"



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
                    CodeString += ", "
            CodeString += ")"

            # check to initialize members with defaultget
            if e.name == g.name:
                CodeString += ": "
                for i in g.introns:
                    CodeString += memprfx+i.name +memsfx+ "("+i.defaultval+")"
                    if i != g.introns[-1]:
                        CodeString += ", "

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
                               type += "* "
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
            CodeString += "#include "+i+";\n\n"
        CodeString += "#include <"+i+".h>\n\n"

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
            CodeString += i.type +memprfx+i.name+memsfx + ";\n"

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
        if e.name == g.name and len(e.introns)==1 and g.name+"& " in e.introns[0].type:
            CodeString += "  ~" + g.name + "();\n"

    #check constructor function definition for private members
    CodeString += "private:\n"
    for i in g.introns:
        if not i.public:
            CodeString += "  "+i.type +memprfx+i.name+memsfx + ";\n"

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
                    i = intron(c, "const char * ", False)
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
                                type += "* "
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

                                create_intron(val, constdef, exon_arg, curr_gene, curr_exon)
                             else:
                                 if ")" in v:
                                     v = v[:v.index(")")]
                                 
                                 # check for a const variable comparison to make private member
                                 create_intron(v, constdef, exon_arg, curr_gene, curr_exon)

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
                            type += "* "
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

    text = file_text
    while True:
        result = re.search('typedef(.*);', text)
        if result:
            val = result.group(1)
            r1 = regex.finditer('\w+$', val)
        
            for r in r1:
                name = r.group(0)
                type =  val[:val.index(name)]
                type = type.lstrip()
                type = type.rstrip()
                g.introns.append(intron(name, type, False))
                break
            text = text[result.start()+len('typedef'+val)+1:]
        else:
            break


def resolve_templates(file_text, g):

    #find name and type, i.e. struct, class
    result = re.search('template[\S\s]*?(?=;)', file_text)
    
    val = result.group(0)
    
    result = re.search('template[\S\s]*?(?=>)>', val)

    if result == None: # error TBD fix!
        print("Discarding rest of file unreadable template discovered....")
        #result = re.search('template[\S\s]*?(?=;)', file_text)
        #file_text = file_text[result.start() + len(result.group(0))+1:]
        file_text = ""
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

  #  print(val)

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
  #  print(e.name)

    for a in args:
        r1 = regex.finditer('\w+$', a)
        for r in r1:
            val = r.group(0)
            # add template intron
            i = intron('<'+val+'>', a[:a.index(val)], False)
            e.introns.append(i)
            break

    # match execution
    result = re.search('template[\S\s]*?(?=;)', file_text)
    val = result.group(0)
    #print(val)
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

        resolve_typedefs(file_text, g)
        resolve_defines(file_text, g)


            #print(file_text)
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
                       # print(file_text)
                        while True:
                            if result2 == None:
                                break

                            file_text = resolve_templates(file_text, g)
                            result2 = re.search('template[\S\s]*?(?=;)', file_text)
                        continue
                    
                    v1 = result.group(0)
                    v1 = v1[:v1.index("(")]

                    if v1 != v:
                       # print("Continuing!")
                        continue

                    if result and result2 and result2.start() < result.start():
                        #print("Continuing!")
                        continue

                    index = result.start()
                    v = v1


                prevfiletext = file_text[:index+len(v)+1]
                file_text = file_text[index:]

                result = re.search(v+'\((.*)\)', file_text)
                if result:
                    args = result.group(1)
                    if ")" in args:
                        args = args[:args.index(")")]
                    args = args.split(",")
                    #val.replace("\n","")

                    if(len(args)==0):
                       
                        result = re.search('(.*)'+v+'\(', prevfiletext)
                        if result:
                            val = result.group(1)
                            #if 'constexpr' in val:
                                #print(dummylines)
                    else:
                        
                        if True:#v=="strlen":
 
                            #print(re.search(v+'\((.*)\)', file_text).start())
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
                                                    #print("NOT FINAL /", val)

                                                   
                                                    #check if final line for define
                                                    txt = file_text[result.start()+len(val):]
                                                    if txt[txt.index('\n')-1] != '\\':
                                                        
                                                        result = re.search(prevregex+'.*\n', file_text)
                                                        val = result.group(0)

                                                        #check to create an exon
                                                        createexon = True
                                                        for e in g.exons:
                                                            if e.name == v.strip():
                                                                createexon = False
                                                                break

                                                        if createexon:
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
                                                #print(v,filename)
                                                 # remove chunk from file text
                                                

                                                #if result == None:
                                                    #print(prevfiletext)
                                                #    print(file_text)
                                                val = result.group(0)
                                                createexon = True
                                                for e in g.exons:
                                                    if e.name == v.strip():
                                                        createexon = False
                                                        break

                                                if createexon:
                                                    e = exon(v.strip())
                                                    e.defenition = val
                                                    g.exons.append(e)

                                                file_text = file_text[result.start()+len(val):]
                                                break
                                elif ' ' in args[0]:
                                    result = re.search('.*?(?='+v+'\()', prevfiletext)
                                    val = result.group(0)
                           
                                    if "#" in val:#skip this line
                                        #print("if! skiping!", file_text)
                                        file_text = file_text[file_text.index('\n'):]
                                        continue
                                     
                                    createexon = True
                                    for e in g.exons:
                                        if e.name == v.strip() and len(args) == len(e.introns):
                                            createexon = False
                                            break

                                    if createexon:
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

                                        #if 'inline ' in e.expression:
                                        #    e.expression = e.expression[e.expression.index("inline ") + len("inline "):]

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
                                                i = intron(val, a[:a.index(val)], False)
                                                e.introns.append(i)
                                                break
                                            #
                                            #print(i.type, i.name)

                                        g.exons.append(e)

                                        #if "__strlen_g" in v:
                                        #print(v,result.group(0))

                                     # try to match funct def
                                                
                                    index = result.start()
                                    #{((?>[^{}]+|(?R))*)}
                                    #{(?:[^{}]+|(?R))*+}

                                    r1 = regex.finditer('{(?:[^{}]+|(?R))*+}', file_text)
                                    result = re.search(v+'\(.*\)\n*\s*\t*?', file_text)
                                    val = result.group(0)
                                    for r in r1:
                                        if r.start() == len(val):
                                           file_text = file_text[r.start() + len(r.group(0)):]
                                        break
                                         
               # print(prevfiletext)
                #if "\n" in prevfiletext:
                #    prevfiletext = prevfiletext[prevfiletext.rindex("\n"):]
                if v+'(' in file_text and file_text.index(v+'(') == 0: 
                    #check to remove funct header
                    result = re.search(v+'\(.*;', file_text)
                    if result:
                        file_text = file_text[result.start() + len(result.group(0)):]
                    else:
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

for m in mutations:
    print(m.name)
    for e in m.exons:
        if e.expression != "":
            val = e.expression + " "+e.name + "("
            for i in e.introns:
                val += i.type + i.name
                if i != e.introns[-1]:
                    val += ", "
            val += ");"
            print(val)



# spawn compiling repos and bring them to life with a thread
index = 0

#while True:
#        
#    #generate individual with base mutation
#    if index < 10:#p.half_limit() != True:
#        fitness = gen_genes(genes, "population/mutant"+str(index), True)
#        t = threading.Thread(name="mutant"+str(index), target=individual, args=("population/mutant"+str(index),fitness, genes,))
#        t.start()
#        index+=1
#
#    if p.extinct():
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