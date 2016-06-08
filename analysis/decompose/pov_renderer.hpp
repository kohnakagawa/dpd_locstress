#pragma once

#include <boost/format.hpp>
#include <string>
#include <sstream>
#include <fstream>

template<class vec3>
std::string make_vector_string(const vec3& vec) {
  std::stringstream ss;
  ss << "<" << vec.x << ", " << vec.y << ", " << vec.z << ">";
  return ss.str();
}

template<class vec3>
struct ObjectParameter {
  std::vector<vec3> color_;
  std::string color_pattern_ = std::string("none");
  vec3 translate_ {0.0, 0.0, 0.0}, scale_ {1.0, 1.0, 1.0}, rotate_ {0.0, 0.0, 0.0};
  
  // texture information
  double transmit_ = 0.0;
  std::string material_ = std::string("none");

  // shape information
  vec3 org_, vec_, cent_;
  double rad_ = 0.5, dist_ = 1.0;
  std::array<vec3, 3> tri_corner_;
  
  // text information
  std::string text_, font_;
  double text_thick_ = 1.0;
};

// color lists
namespace PovColor {
  template<class vec3>
  constexpr vec3 black() {
    return {0.0, 0.0, 0.0};
  }

  template<class vec3>
  constexpr vec3 white() {
    return {1.0, 1.0, 1.0};
  }

  template<class vec3>
  constexpr vec3 red() {
    return {1.0, 0.0, 0.0};
  }

  template<class vec3>
  constexpr vec3 green() {
    return {0.0, 1.0, 0.0};
  }

  template<class vec3>
  constexpr vec3 blue() {
    return {0.0, 0.0, 1.0};
  }

  template<class vec3>
  constexpr vec3 grey() {
    return {0.5, 0.5, 0.5};
  }

  template<class vec3>
  constexpr vec3 yellow() {
    return {1.0, 1.0, 0.0};
  }

  template<class vec3>
  constexpr vec3 cyan() {
    return {0.0, 1.0, 1.0};
  }

  template<class vec3>
  constexpr vec3 magenta() {
    return {1.0, 0.0, 1.0};
  }

  template<class vec3>
  constexpr vec3 orange_yellow() {
    return {1.0, 0.65, 0.0};
  }

  template<class vec3>
  constexpr vec3 yellow_green() {
    return {0.25, 1.0, 0.0};
  }

  template<class vec3>
  constexpr vec3 red_wine() {
    return {1.0, 0.0, 0.25};
  }

  template<class vec3>
  constexpr vec3 dark_purple() {
    return {0.53, 0.12, 0.47};
  }
}

// interface class
template<class vec3>
class DrawnObject {
protected:
  std::vector<vec3> color_;
  std::string color_pattern_;
  std::string material_;
  vec3 translate_, scale_, rotate_;
  double transmit_;
  std::string obj_modif_;

  void SetObjectModifiers() {
    std::stringstream ss;
    if (material_ != "none") {
      ss << boost::format("material {%s} scale %s translate %s rotate %s transmit %f")
	% material_ % make_vector_string(scale_) % make_vector_string(translate_) % make_vector_string(rotate_) % transmit_;
    } else {
      if (color_pattern_ == "none") {
	ss << boost::format("texture {pigment {color rgb %s} finish {phong 1}} scale %s translate %s rotate %s transmit %f")
	  % make_vector_string(color_[0]) % make_vector_string(scale_) % make_vector_string(translate_) % make_vector_string(rotate_) % transmit_;      
      } else {
	// ASSUME color_.size() == 2
	ss << boost::format("texture {pigment {%s color rgb %s color rgb %s} finish {phong 1}} scale %s translate %s rotate %s transmit %f")
	  % color_pattern_ % make_vector_string(color_[0]) % make_vector_string(color_[1])
	  % make_vector_string(scale_) % make_vector_string(translate_) % make_vector_string(rotate_) % transmit_;
      }
    }
    obj_modif_ = ss.str();
  }
  
public:
  DrawnObject(const ObjectParameter<vec3>& param) {
    this->color_	 = param.color_;
    this->color_pattern_ = param.color_pattern_;
    this->material_      = param.material_;
    this->translate_	 = param.translate_;
    this->rotate_	 = param.rotate_;
    this->scale_	 = param.scale_;
    this->transmit_      = param.transmit_;
    SetObjectModifiers();
  }
  virtual void WriteToOstream(std::ostream& ost) = 0;
  virtual ~DrawnObject() = default;
};

template<class vec3>
class VectorPov : public DrawnObject<vec3> {
  vec3 org_, vec_;
  double rad_;
  
  void WriteVectorPovMacro(std::ostream& ost) {
    ost << "#macro Vector (start_p, end_p, radius_cyl)\n";
    ost << "union {\n";
    ost << "  cylinder {start_p, end_p - (vnormalize(end_p - start_p) * 4.25 * radius_cyl), radius_cyl}\n";
    ost << "  cone {end_p - (vnormalize(end_p - start_p) * 5.0 * radius_cyl), 3 * radius_cyl, end_p, 0.0}\n";
    ost << "} // end of union\n";
    ost << "#end // end of macro\n\n";
  }
  
public:
  VectorPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->org_ = param.org_;
    this->vec_ = param.vec_;
    this->rad_ = param.rad_;
  }
  ~VectorPov() = default;

  void WriteToOstream(std::ostream& ost) override {
    static bool is_first_call = true;
    if (is_first_call) {
      WriteVectorPovMacro(ost);
      is_first_call = false;
    }
    
    const auto end_p = org_ + vec_;
    ost << boost::format("object {Vector(%s, %s, %f) %s}\n\n")
      % make_vector_string(org_) % make_vector_string(end_p) % rad_ % this->obj_modif_;
  }
};

template<class vec3>
class CylinderPov : public DrawnObject<vec3> {
  vec3 org_, vec_;
  float rad_;
public:
  CylinderPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->org_ = param.org_;
    this->vec_ = param.vec_;
    this->rad_ = param.rad_;
  }
  
  void WriteToOstream(std::ostream& ost) override {
    const auto end_pos = org_ + vec_;
    ost << boost::format("cylinder {%s, %s, %f %s}\n\n")
      % make_vector_string(org_) % make_vector_string(end_pos) % rad_ % this->obj_modif_;
  }
};

template<class vec3>
class TriSurfPov : public DrawnObject<vec3> {
  std::array<vec3, 3> tri_corner_;
public:
  TriSurfPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->tri_corner_ = param.tri_corner_;
  }
  
  void WriteToOstream(std::ostream& ost) override {
    ost << boost::format("triangle {%s, %s, %s %s}\n\n")
      % make_vector_string(tri_corner_[0]) % make_vector_string(tri_corner_[1]) % make_vector_string(tri_corner_[2]) % this->obj_modif_;
  }
};

template<class vec3>
class PlanePov : public DrawnObject<vec3> {
  vec3 vec_;
  double dist_;
public:
  PlanePov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->vec_  = param.vec_;
    this->dist_ = param.dist_;
  }

  void WriteToOstream(std::ostream& ost) override {
    ost << boost::format("plane {%s, %f %s}\n\n")
      % make_vector_string(vec_) % dist_ % this->obj_modif_;
  }
};

template<class vec3>
class SpherePov : public DrawnObject<vec3> {
  vec3 cent_;
  double rad_;
public:
  SpherePov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->cent_ = param.cent_;
    this->rad_  = param.rad_;
  }

  void WriteToOstream(std::ostream& ost) override {
    ost << boost::format("sphere {%s, %f %s}\n\n")
      % make_vector_string(cent_) % rad_ % this->obj_modif_;
  }
};

template<class vec3>
class TextPov : public DrawnObject<vec3> {
  std::string text_, font_;
  double thick_ = 1.0;
public:
  TextPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->text_ = param.text_;
    this->font_ = param.font_;
    this->thick_ = param.text_thick_;
  }

  void WriteToOstream(std::ostream& ost) override {
    std::stringstream ss;
    ss << boost::format("texture {pigment {color Black}} scale %s translate %s rotate %s transmit %f")
      % make_vector_string(this->scale_) % make_vector_string(this->translate_) % make_vector_string(this->rotate_) % this->transmit_;
    this->obj_modif_ = ss.str();
    ost << boost::format("object{text {ttf \"%s\", \"%s\", %f, 0 %s}}\n\n")
      % font_ % text_ % thick_ % this->obj_modif_;
  }
};

#define SET_COORD(param)			\
  do {						\
    param.translate_	= translate;		\
    param.scale_	= scale;		\
    param.rotate_	= rotate;		\
  } while (false)

template<class vec3>
class ObjectFactory {
public:
  static DrawnObject<vec3>* CreateVector(const vec3& org,
					 const vec3& vec,
					 const double rad,
					 const vec3& color,
					 const vec3& translate,
					 const vec3& scale,
					 const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.org_		= org;
    param.vec_		= vec;
    param.rad_          = rad;
    param.color_.push_back(color);
    SET_COORD(param);
    return new VectorPov<vec3>(param);
  }

  // Do not translate or rotate
  static DrawnObject<vec3>* CreateVector(const vec3& org,
					 const vec3& vec,
					 const double rad,
					 const vec3& color) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreateVector(org, vec, rad, color, translate, scale, rotate);
  }
  
  static DrawnObject<vec3>* CreateCylinder(const vec3& org,
					   const vec3& vec,
					   const double rad,
					   const vec3& color,
					   const vec3& translate,
					   const vec3& scale,
					   const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.org_		= org;
    param.vec_		= vec;
    param.rad_		= rad;
    param.color_.push_back(color);
    SET_COORD(param);
    return new CylinderPov<vec3>(param);
  }
  
  // Do not translate or rotate
  static DrawnObject<vec3>* CreateCylinder(const vec3& org,
					   const vec3& vec,
					   const vec3& color,
					   const double rad) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreateCylinder(org, vec, rad, color, translate, scale, rotate);
  }
  
  // set material property
  static DrawnObject<vec3>* CreateCylinder(const vec3& org,
					   const vec3& vec,
					   const vec3& color,
					   const std::string material,
					   const double rad) {
    ObjectParameter<vec3> param;
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    param.color_.push_back(color);
    SET_COORD(param);
    param.material_ = material;
    return CreateCylinder(org, vec, rad, color, translate, scale, rotate);
  }

  
  static DrawnObject<vec3>* CreateTriSurf(const std::array<vec3, 3>& tri_corner,
					  const vec3& color,
					  const vec3& translate,
					  const vec3& scale,
					  const vec3& rotate,
					  const double transmit) {
    ObjectParameter<vec3> param;
    param.tri_corner_	= tri_corner;
    param.color_.push_back(color);
    SET_COORD(param);
    param.transmit_ = transmit;
    param.material_ = "M_Glass";
    return new TriSurfPov<vec3>(param);
  }

  // Do not translate or rotate
  static DrawnObject<vec3>* CreateTriSurf(const std::array<vec3, 3>& tri_corner,
					  const vec3& color,
					  const double transmit) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreateTriSurf(tri_corner, color, translate, scale, rotate, transmit);
  }

  static DrawnObject<vec3>* CreatePlane(const vec3& vec,
					const double dist,
					const std::vector<vec3>& color,
					const std::string pattern,
					const vec3& translate,
					const vec3& scale,
					const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.vec_ = vec;
    param.dist_ = dist;
    param.color_ = color;
    param.color_pattern_ = pattern;
    SET_COORD(param);
    return new PlanePov<vec3>(param);
  }
  
  static DrawnObject<vec3>* CreatePlane(const vec3& vec,
					const double dist,
					const std::vector<vec3>& color,
					const std::string pattern) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreatePlane(vec, dist, color, pattern, translate, scale, rotate);
  }

  static DrawnObject<vec3>* CreateSphere(const vec3& cent,
					 const double rad,
					 const vec3& color,
					 const vec3& translate,
					 const vec3& scale,
					 const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.cent_ = cent;
    param.rad_  = rad;
    param.color_.push_back(color);
    SET_COORD(param);
    return new SpherePov<vec3>(param);
  }

  // Do not translate or rotate
  static DrawnObject<vec3>* CreateSphere(const vec3& cent,
					 const double rad,
					 const vec3& color) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreateSphere(cent, rad, color, translate, scale, rotate);
  }
  
  static DrawnObject<vec3>* CreateText(const std::string text,
				       const std::string font,
				       const double text_thick,
				       const vec3& color,
				       const vec3& translate,
				       const vec3& scale,
				       const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.text_	 = text;
    param.font_	 = font;
    param.text_thick_ = text_thick;
    param.color_.push_back(color);
    SET_COORD(param);
    return new TextPov<vec3>(param);
  }
  
  // Do not translate or rotate
  static DrawnObject<vec3>* CreateText(const std::string text,
				       const std::string font,
				       const double text_thick,
				       const vec3& color) {
    const vec3 translate {0.0, 0.0, 0.0};
    const vec3 scale {1.0, 1.0, 1.0};
    const vec3 rotate {0.0, 0.0, 0.0};
    return CreateText(text, font, text, color, translate, scale, rotate);
  }  
};

#undef SET_COORD

template<class vec3>
class PovRenderer {
  std::ofstream fout;

  std::vector<std::unique_ptr<DrawnObject<vec3>>> ptr_objects;
  
  std::vector<std::string> include_file_list;
  std::string camera_settings;
  std::string background_settings;
  std::string light_settings;

  void AddPreprocessToScript() {
    for (const auto& file : include_file_list)
      fout << "#include \"" + file + "\"\n";
  }

public:
  PovRenderer(const char* fname) {
    fout.open(fname);
  }

  void SetCamera(const vec3& camera_loc,
		 const vec3& camera_look_at) {
    const auto location = "location " + make_vector_string(camera_loc) + " ";
    const auto look_at  = "look_at "  + make_vector_string(camera_look_at);
    camera_settings = "camera {" + location + look_at + "}\n\n";
  }

  void SetBackGround() {
    // NOTE: always white
    background_settings = "background {color White}\n\n";
  }

  void SetLight(const vec3& light_at) {
    light_settings = "light_source {" + make_vector_string(light_at) + " color White}\n\n";
  }
  
  void AppendObject(DrawnObject<vec3>* ptr) {
    ptr_objects.push_back(std::unique_ptr<DrawnObject<vec3>>(ptr));
  }

  void AppendIncludeFile(const std::string file) {
    include_file_list.push_back(file);
  }

  void AppendDefaultIncludeFiles() {
    AppendIncludeFile("colors.inc");
    AppendIncludeFile("shapes.inc");
    AppendIncludeFile("textures.inc");
    AppendIncludeFile("glass.inc");
    AppendIncludeFile("metals.inc");
  }

  void WriteRenderScript() {
    // set include files
    AddPreprocessToScript();
    // set background
    fout << background_settings;
    // set light
    fout << light_settings;
    // set camera
    fout << camera_settings;
    // draw objects
    for (auto& ptr_obj : ptr_objects) {
      ptr_obj->WriteToOstream(fout);
    }
  }
};

