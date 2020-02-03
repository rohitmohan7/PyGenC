import sys
import subprocess
import os
import re
import shutil

memprfx = ""
memsfx = "_"
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

def create_cpp(g, path):

    print("generating cpp... ", g.name +(".cpp"))
    CodeString = "\n/**************** Auto Generated File **********************/\n\n"

    #include header
    CodeString += "#include \"../includes/" +g.name+".h\"\n\n"
    
    #define Codons?

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
            for s in e.set:
                # check if constructor
                if e.name == g.name:
                    #check for another function that may be setting this 
                    for e2 in g.exons:
                        set = False
                        if e2.name != g.name:
                            for s2 in e2.set:
                                if s2.name == s.name:
                                    CodeString += "   "+e2.name+"("+s.name+");\n"
                                    set = True
                                    break
                            if set:
                                break
                else:
                     if "const char *" in s.type:
                        
                        CodeString += "  if ("+s.name+")\n  {\n"
                        CodeString += "    const size_t len = strlen("+s.name + ");\n"
                        CodeString += "    char* const clone = new char[len + 1];\n"
                        CodeString += "    memcpy(clone, "+s.name+", len + 1);\n"
                        CodeString += "\n    if ("+memprfx + s.name + memsfx +" != NULL)\n    {\n"
                        CodeString += "      delete[] "+memprfx + s.name + memsfx + ";\n    }\n\n"
                        CodeString += "    "+memprfx + s.name + memsfx +" = clone;\n" 
                        CodeString += "  }\n  else\n  {\n"
                        CodeString += "    delete[] "+memprfx + s.name + memsfx + ";\n"
                        CodeString += "    "+memprfx + s.name + memsfx + " = NULL;\n  }\n"

                        #add to includes 
                        if not "string" in g.includes:
                            g.includes.append("string")

            if e.get:

                for i in g.introns:
                    if i.name == e.get.name:
                        # check return type to return member here within all C++ rules and make a sensible conversion
                        if e.expression in i.type:
                            CodeString += "   return "+memprfx+e.get.name+memsfx+";\n"
                        elif "int" in e.expression and "const char *" in i.type: # send length
                            CodeString += "   if(!"+memprfx+i.name+memsfx+")\n    return "+e.defaultget+";\n"
                            CodeString += "   return strlen("+memprfx+i.name+memsfx+");\n"
                            if not "string" in g.includes:
                                g.includes.append("string")
                        break

            # define copy constructor
            if(e.name == g.name and len(e.introns) == 1 and g.name+"& " in e.introns[0].type):            
                # find a method to copy from rhs all introns
                for i in g.introns:
                    method = False
                    for e2 in g.exons:
                        if e2.name != g.name:
                            for s in e2.set:
                                if s.name == i.name:
                                    CodeString += "  "+e2.name+"("+e.introns[0].name+"."+memprfx + i.name + memsfx+");\n"
                                    method = True
                                    break
                        if method:
                            break
                    if not method:
                        CodeString += "  " +memprfx + i.name + memsfx+ " = "+e.introns[0].name +"."+memprfx + i.name + memsfx+";\n"
                CodeString += "}\n\n"

                # define deconstructor
                CodeString += g.name + "::~"+g.name+"()\n{\n"
                for i in g.introns:
                    if "*" in i.type:
                        CodeString += "  delete[] "+memprfx+i.name+memsfx+";\n"

            CodeString += "}\n\n"

    cfile = open(path+"/src/"+ g.name + ".cpp", 'w')
    cfile.write(CodeString)
    cfile.close()
    return

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
            CodeString += i.type +memprfx+i.name+memsfx + ";\n";

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

        CodeString += ")";
        # constructor member intialization in cpp
        CodeString += ";\n";

        #check for deconstructor for copy constructor, need to check if virtualization is needed!
        if e.name == g.name and len(e.introns)==1 and g.name+"& " in e.introns[0].type:
            CodeString += "  ~" + g.name + "();\n"

    #check constructor function definition for private members
    CodeString += "private:\n"
    for i in g.introns:
        if not i.public:
            CodeString += "  "+i.type +memprfx+i.name+memsfx + ";\n";

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

for g in genes:

    # construct base repo dir
    construct_repo_dir("base_repo")

    # create definitions
    create_cpp(g, "./base_repo")
    create_header(g, "./base_repo")
    create_test(g, "./base_repo")
        
if not os.path.exists(os.path.dirname("./base_repo/")):
    print("Failed to create base_repo, check Testdef")
    exit()
print("compiling base_repo")
try:
    create_testmain(genes, "./base_repo")
    create_cmakelists(genes, "./base_repo")
    
    cmakeCmd = ["cmake", "base_repo/CMakeLists.txt"]
    retCode = subprocess.check_call(cmakeCmd)
    
    cmakeCmd = ["make", "-C", "base_repo"]
    retCode = subprocess.check_call(cmakeCmd)

    print("compiled base_repo........executing tests")

except subprocess.CalledProcessError:
    print("failed to compile base_repo")
    exit()
except OSError:
    print("failed to compile base_repo")
    exit()
    # create defenition for constructor 
cmakeCmd = ["./base_repo/executeTests"]
retCode = subprocess.check_call(cmakeCmd)
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