
include conf/mara.conf

MARA_E  = bin/mara
MARA_A  = lib/libmara.a
LUA_A   = lib/liblua.a
LUNUM_A = lib/liblunum.a
GLFW_A  = lib/libglfw.a

ifeq ($(USE_GLFW), False)
GLFW_A =
GL_L   =
endif

LUA_VERSION  = lua-5.2.0
GLFW_VERSION = glfw-2.7.2
LUNUM_BRANCH = master
LUNUM_GIT    = git@github.com:jzrake/lunum.git
SERVE_SRC    = http://zrake.webfactional.com/media



# ------------------------------------------------------------------------------
# Mara executable
# ------------------------------------------------------------------------------
$(MARA_E) : $(LUA_A) $(LUNUM_A) $(GLFW_A) $(MARA_A)
	$(CXX) -o $@ -Llib -lmara -llunum $(LUA_A) $(HDF5_L) $(FFTW_L) $(GL_L) $(CLIBS)



# ------------------------------------------------------------------------------
# Library build aliases
# ------------------------------------------------------------------------------
mara   :  $(MARA_A)
lua    :  $(LUA_A)
lunum  :  $(LUNUM_A)
glfw   :  $(GLFW_A)


$(MARA_A) : FORCE
	@make -C src

$(LUNUM_A) :
	@echo "Downloading and installing Lunum"
	@mkdir -p bin lib include
	git clone -b $(LUNUM_BRANCH) $(LUNUM_GIT)
	make -C lunum install CC=$(CC) LUA_HOME=$(PWD) INSTALL_TOP=$(PWD)
	rm -rf lunum

$(LUA_A) : 
	@echo "Downloading and installing Lua..."
	@mkdir -p bin lib include
	wget $(SERVE_SRC)/$(LUA_VERSION).tar.gz
	tar xvf $(LUA_VERSION).tar.gz
	make -C $(LUA_VERSION) $(ARCH_LUA) install CC=$(CC) INSTALL_TOP=$(PWD)
	cp src/*.lua lib/lua/5.2
	rm -rf share man $(LUA_VERSION) $(LUA_VERSION).tar.gz

$(GLFW_A) :
	@echo "Downloading and installing glfw..."
	@mkdir -p bin lib include
	wget $(SERVE_SRC)/$(GLFW_VERSION).tar.gz
	tar xvf $(GLFW_VERSION).tar.gz
	make -C $(GLFW_VERSION) $(ARCH_GLFW) CC=$(CC)
	cp -r $(GLFW_VERSION)/lib/$(ARCH_GLFW)/libglfw.a lib
	cp -r $(GLFW_VERSION)/include/* include
	rm -r $(GLFW_VERSION) $(GLFW_VERSION).tar.gz


# Clean up Mara src and executable
clean : 
	@make -C src clean
	@rm -f bin/mara

# Clean up dependincies as well: lib, bin, include directories
realclean : clean
	@rm -rf bin lib include

FORCE :
