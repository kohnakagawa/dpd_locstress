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

// object parameter class
template<class vec3>
class ObjectParameter {
public:
  vec3 color_ {0.6, 1.0, 0.0};
  // TODO: add object texture
  vec3 translate_ {0.0, 0.0, 0.0}, scale_ {1.0, 1.0, 1.0}, rotate_ {0.0, 0.0, 0.0};

  // shape information
  vec3 org_, vec_;
  double rad_ = 0.5;
  std::array<vec3, 3> tri_corner_;
  
  // text information
  std::string text_, font_;
  double text_thick_ = 1.0;

  // texture information
  
};

// interface class
template<class vec3>
class DrawnObject {
protected:
  vec3 color_;
  vec3 translate_, scale_, rotate_;
  std::string obj_modif_;

  void SetObjectModifiers() {
    std::stringstream ss;
    ss << boost::format("pigment {color rgb %s} translate %s scale %s rotate %s")
      % make_vector_string(color_) % make_vector_string(translate_) % make_vector_string(scale_) % make_vector_string(rotate_);
    obj_modif_ = ss.str();
  }
  
public:
  DrawnObject(const ObjectParameter<vec3>& param) {
    this->color_	= param.color_;
    this->translate_	= param.translate_;
    this->rotate_	= param.rotate_;
    this->scale_	= param.scale_;
    SetObjectModifiers();
  }
  virtual void WriteToOstream(std::ostream& ost) = 0;
  virtual ~DrawnObject() = default;
};

template<class vec3>
class VectorPov : public DrawnObject<vec3> {
  vec3 org_, vec_;
  const float cyl_rad = 0.2;
  
  void WriteVectorPovMacro(std::ostream& ost) {
    ost << "#macro Vector (start_p, end_p, radius_cyl)\n";
    ost << "union {\n";
    ost << "  cylinder {start_p, end_p - (vnormalize(end_p - start_p) * 9.5 * radius_cyl), radius_cyl}\n";
    ost << "  cone {end_p - (vnormalize(end_p - start_p) * 10.0 * radius_cyl), 3 * radius_cyl, end_p, 0.0}\n";
    ost << "} // end of union\n";
    ost << "#end // end of macro\n\n";
  }
  
public:
  VectorPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->org_ = param.org_;
    this->vec_ = param.vec_;
  }
  ~VectorPov() = default;

  void WriteToOstream(std::ostream& ost) override {
    static bool is_first_call = true;
    if (is_first_call) {
      WriteVectorPovMacro(ost);
      is_first_call = false;
    }
    
    ost << boost::format("object {Vector(%s, %s, %f) %s}\n\n")
      % make_vector_string(org_) % make_vector_string(vec_) % cyl_rad % this->obj_modif_;
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
class TextPov : public DrawnObject<vec3> {
  std::string text_, font_;
  float thick_ = 1.0;
public:
  TextPov(const ObjectParameter<vec3>& param) : DrawnObject<vec3>(param) {
    this->text_ = param.text_;
    this->font_ = param.font_;
    this->thick_ = param.text_thick_;
  }

  void WriteToOstream(std::ostream& ost) override {
    ost << boost::format("text {ttf \"%s\", \"%s\", %f, 0 %s}\n\n")
      % text_ % font_ % thick_ % this->obj_modif_;
  }
};

#define SET_COLOR_AND_COORD(param)		\
  do {						\
    param.color_	= color;		\
    param.translate_	= translate;		\
    param.scale_	= scale;		\
    param.rotate_	= rotate;		\
  } while (false)

template<class vec3>
class ObjectFactory {
public:
  static DrawnObject<vec3>* CreateVector(const vec3& org,
					 const vec3& vec,
					 const vec3& color,
					 const vec3& translate,
					 const vec3& scale,
					 const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.org_		= org;
    param.vec_		= vec;
    SET_COLOR_AND_COORD(param);
    return new VectorPov<vec3>(param);
  }

  // Do not translate or rotate
  static DrawnObject<vec3>* CreateVector(const vec3& org,
					 const vec3& vec,
					 const vec3& color) {
    
    ObjectParameter<vec3> param;
    param.org_		= org;
    param.vec_		= vec;
    param.color_        = color;
    return new VectorPov<vec3>(param);
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
    SET_COLOR_AND_COORD(param);
    return new CylinderPov<vec3>(param);
  }
  
  // Do not translate or rotate
  static DrawnObject<vec3>* CreateCylinder(const vec3& org,
					   const vec3& vec,
					   const vec3& color) {
    ObjectParameter<vec3> param;
    param.org_		= org;
    param.vec_		= vec;
    param.color_        = color;
    return new CylinderPov<vec3>(param);
  }
  
  static DrawnObject<vec3>* CreateTriSurf(const std::array<vec3, 3>& tri_corner,
					  const vec3& color,
					  const vec3& translate,
					  const vec3& scale,
					  const vec3& rotate) {
    ObjectParameter<vec3> param;
    param.tri_corner_	= tri_corner;
    SET_COLOR_AND_COORD(param);
    return new TriSurfPov<vec3>(param);
  }

  // Do not translate or rotate
  static DrawnObject<vec3>* CreateTriSurf(const std::array<vec3, 3>& tri_corner,
					  const vec3& color) {
    ObjectParameter<vec3> param;
    param.tri_corner_	= tri_corner;
    param.color_	= color;
    return new TriSurfPov<vec3>(param);    
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
    SET_COLOR_AND_COORD(param);
    return new TextPov<vec3>(param);
  }
  
  // Do not translate or rotate
  static DrawnObject<vec3>* CreateText(const std::string text,
				       const std::string font,
				       const double text_thick,
				       const vec3& color) {
    ObjectParameter<vec3> param;
    param.text_	 = text;
    param.font_	 = font;
    param.text_thick_ = text_thick;
    return new TextPov<vec3>(param);
  }  
};

#undef SET_COLOR_AND_COORD

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

  void SetCameraSettings(const vec3& camera_loc,
			 const vec3& camera_look_at) {
    const auto location = "location " + make_vector_string(camera_loc)     + " ";
    const auto look_at  = "look_at "  + make_vector_string(camera_look_at);
    camera_settings = "camera {" + location + look_at + "}\n\n";
  }

  void SetBackGround() {
    background_settings = "background {color White}\n\n";
  }

  void SetLight(const vec3& light_at) {
    light_settings = "light_source {" + make_vector_string(light_at) + " color White}\n\n";
  }
  
  void AppendObject(DrawnObject<vec3>* ptr) {
    ptr_objects.push_back(std::unique_ptr<DrawnObject<vec3>>(ptr));
  }

  void AppendIncludeFiles(const std::string file) {
    include_file_list.push_back(file);
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
    for (size_t i = 0; i < ptr_objects.size(); i++)
      ptr_objects[i]->WriteToOstream(fout);
  }
};
